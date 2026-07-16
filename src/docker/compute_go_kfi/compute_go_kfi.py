#!/usr/bin/env python3
"""
GO Knock-out Feasibility Index (GO-KFI) 計算スクリプト
各GO termに対してノックアウトの実行可能性を評価
"""

import sys
import pandas as pd
import numpy as np
import argparse


def _parse_guide_target_set(guide_targets_str):
    """Return unique, non-empty guide targets from a serialized field."""
    if pd.isna(guide_targets_str):
        return set()

    guide_field = str(guide_targets_str).strip()
    if not guide_field or guide_field.lower() == 'nan':
        return set()

    return {
        guide.strip()
        for guide in guide_field.split(',')
        if guide.strip()
    }


def parse_guide_targets(guide_targets_str):
    """Return the number of unique guide targets in a serialized field."""
    return len(_parse_guide_target_set(guide_targets_str))


def _require_valid_unique_ids(df):
    """Require every annotation row to have a usable gene-level identifier."""
    if 'unique_id' not in df.columns:
        raise ValueError("Required column 'unique_id' was not found.")

    normalized_ids = df['unique_id'].astype('string').str.strip()
    invalid_ids = (
        normalized_ids.isna()
        | normalized_ids.eq('')
        | normalized_ids.str.lower().eq('nan')
    )

    invalid_count = int(invalid_ids.sum())
    if invalid_count > 0:
        raise ValueError(
            f"Found {invalid_count} annotation rows without a valid unique_id."
        )


def _parse_guide_count(value, unique_id):
    """Return one validated non-negative integer guide count."""
    if pd.isna(value) or isinstance(value, (bool, np.bool_)):
        raise ValueError(
            f"guide_count for {unique_id!r} must be a non-negative integer."
        )

    if isinstance(value, (int, np.integer)):
        guide_count = int(value)
    elif isinstance(value, (float, np.floating)):
        if not np.isfinite(value) or not float(value).is_integer():
            raise ValueError(
                f"guide_count for {unique_id!r} must be a non-negative integer."
            )
        guide_count = int(value)
    else:
        normalized = str(value).strip()
        if not normalized.isdigit():
            raise ValueError(
                f"guide_count for {unique_id!r} must be a non-negative integer."
            )
        guide_count = int(normalized)

    if guide_count < 0:
        raise ValueError(
            f"guide_count for {unique_id!r} must be a non-negative integer."
        )
    return guide_count


def _validated_guide_counts_by_unique_id(df):
    """Validate canonical guide counts and reject per-gene conflicts."""
    guide_counts = {}
    for unique_id_value, guide_count_value in zip(
        df['unique_id'], df['guide_count']
    ):
        unique_id = str(unique_id_value).strip()
        guide_count = _parse_guide_count(guide_count_value, unique_id)
        previous_count = guide_counts.get(unique_id)
        if previous_count is not None and previous_count != guide_count:
            raise ValueError(
                f"Conflicting guide_count values for {unique_id!r}: "
                f"{previous_count} and {guide_count}."
            )
        guide_counts[unique_id] = guide_count
    return guide_counts


def _merge_gene_guides(gene_data, guide_targets_str):
    """Merge guide targets while avoiding large sets for identical repeated rows."""
    if pd.isna(guide_targets_str):
        return

    guide_field = str(guide_targets_str).strip()
    if not guide_field or guide_field.lower() == 'nan':
        return

    if gene_data['guide_set'] is not None:
        gene_data['guide_set'].update(_parse_guide_target_set(guide_field))
    elif gene_data['guide_field'] is None:
        gene_data['guide_field'] = guide_field
    elif gene_data['guide_field'] != guide_field:
        gene_data['guide_set'] = _parse_guide_target_set(
            gene_data['guide_field']
        )
        gene_data['guide_set'].update(_parse_guide_target_set(guide_field))
        gene_data['guide_field'] = None


def _get_gene_guide_count(gene_data):
    """Return a gene's unique guide count after row-level aggregation."""
    if gene_data['guide_count'] is not None:
        return gene_data['guide_count']
    if gene_data['guide_set'] is not None:
        return len(gene_data['guide_set'])
    return parse_guide_targets(gene_data['guide_field'])


def calculate_go_kfi(df):
    """Calculate GO-KFI using unique_id as the gene-level counting unit."""
    _require_valid_unique_ids(df)
    canonical_guide_counts = (
        _validated_guide_counts_by_unique_id(df)
        if 'guide_count' in df.columns
        else None
    )

    genes = {}
    go_order = []
    seen_go_terms = set()

    # Aggregate annotation rows into genes before calculating GO statistics.
    for _, row in df.iterrows():
        unique_id = str(row['unique_id']).strip()
        gene_data = genes.setdefault(unique_id, {
            'go_terms': set(),
            'guide_count': (
                canonical_guide_counts[unique_id]
                if canonical_guide_counts is not None
                else None
            ),
            'guide_field': None,
            'guide_set': None,
        })

        if canonical_guide_counts is None:
            _merge_gene_guides(gene_data, row.get('guide_targets', ''))

        go_terms_str = row.get('GO_ALL', '')
        if pd.isna(go_terms_str) or not str(go_terms_str).strip():
            continue

        for go_term in str(go_terms_str).split(';'):
            go_term = go_term.strip()
            if not go_term:
                continue

            gene_data['go_terms'].add(go_term)
            if go_term not in seen_go_terms:
                seen_go_terms.add(go_term)
                go_order.append(go_term)

    go_stats = {
        go_term: {
            'total_genes': 0,
            'genes_with_guides': 0,
            'total_guides': 0,
        }
        for go_term in go_order
    }

    # Add each gene to each associated GO term exactly once.
    for gene_data in genes.values():
        guide_count = _get_gene_guide_count(gene_data)
        for go_term in gene_data['go_terms']:
            stats = go_stats[go_term]
            stats['total_genes'] += 1
            stats['total_guides'] += guide_count
            if guide_count > 0:
                stats['genes_with_guides'] += 1

    go_kfi_results = []

    for go_term in go_order:
        stats = go_stats[go_term]
        total_genes = stats['total_genes']
        genes_with_guides = stats['genes_with_guides']
        total_guides = stats['total_guides']

        if total_genes > 0:
            kfi_score = genes_with_guides / total_genes
        else:
            kfi_score = 0.0

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

    output_columns = [
        'GO_ID',
        'total_genes',
        'genes_with_guides',
        'total_guides',
        'kfi_score',
        'targetable_type',
    ]
    return pd.DataFrame(go_kfi_results, columns=output_columns)


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


def validate_input_data(df):
    """入力データの検証"""
    print("Validating input data...")

    try:
        _require_valid_unique_ids(df)
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return False
    
    # GO列の確認
    go_columns = [col for col in df.columns if col.startswith('GO_')]
    print(f"Found GO columns: {go_columns}")
    
    if not go_columns:
        print("ERROR: No GO columns found in input data!", file=sys.stderr)
        print("Expected columns starting with 'GO_'", file=sys.stderr)
        return False
    
    # GO列にデータが存在するかチェック
    total_go_entries = 0
    for go_col in go_columns:
        non_empty = df[go_col].notna() & (df[go_col] != '') & (df[go_col] != 'nan')
        count = non_empty.sum()
        total_go_entries += count
        print(f"  {go_col}: {count} non-empty entries ({count/len(df)*100:.1f}%)")
    
    if total_go_entries == 0:
        print("ERROR: No GO annotation data found!", file=sys.stderr)
        print("All GO columns are empty. This indicates GO annotation retrieval failed.", file=sys.stderr)
        print("Possible causes:", file=sys.stderr)
        print("1. UniProt API connection failure", file=sys.stderr)
        print("2. Invalid UniProt IDs in subject_id column", file=sys.stderr)
        print("3. Network connectivity issues", file=sys.stderr)
        return False
    
    if 'guide_count' in df.columns:
        try:
            canonical_guide_counts = _validated_guide_counts_by_unique_id(df)
        except ValueError as error:
            print(f"ERROR: {error}", file=sys.stderr)
            return False
        genes_with_guides = sum(
            guide_count > 0 for guide_count in canonical_guide_counts.values()
        )
        total_genes = len(canonical_guide_counts)
        percentage = (
            genes_with_guides / total_genes * 100 if total_genes else 0.0
        )
        print(
            "guide_count: "
            f"{genes_with_guides} genes with guides ({percentage:.1f}%)"
        )
        if genes_with_guides == 0:
            print("WARNING: No guide targets found!", file=sys.stderr)
    else:
        # Preserve support for contig-level legacy tables.
        if 'guide_targets' not in df.columns:
            print(
                "ERROR: guide_count or guide_targets column not found!",
                file=sys.stderr,
            )
            return False

        guide_non_empty = (
            df['guide_targets'].notna()
            & (df['guide_targets'] != '')
            & (df['guide_targets'] != 'nan')
        )
        guide_count = guide_non_empty.sum()
        print(
            f"guide_targets: {guide_count} non-empty entries "
            f"({guide_count/len(df)*100:.1f}%)"
        )

        if guide_count == 0:
            print("WARNING: No guide targets found!", file=sys.stderr)
    
    print(f"Data validation completed. Total GO entries: {total_go_entries}")
    return True

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
        print(f"Successfully loaded {len(df)} rows with {len(df.columns)} columns")
        print(f"Columns: {list(df.columns)}")
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # 入力データの検証
    if not validate_input_data(df):
        print("Data validation failed. Cannot proceed with GO-KFI calculation.", file=sys.stderr)
        
        # 空の結果ファイルを作成
        empty_df = pd.DataFrame(columns=['GO_ID', 'total_genes', 'genes_with_guides', 'total_guides', 'kfi_score', 'targetable_type'])
        empty_df.to_csv(args.output, index=False)
        
        # エラー情報を含む統計ファイルを作成
        error_stats = """[ERROR] GO-KFI calculation failed because input validation failed.

Root cause: Required GO annotations, guide targets, or unique gene identifiers are missing or invalid.

Recommendations:
1. Check network connectivity to UniProt REST API
2. Verify UniProt IDs in subject_id column are valid
3. Verify every annotation row has a non-empty unique_id
4. Consider using a local GO annotation database
5. Re-run the go_analysis step with valid input data

Statistics: 0 GO terms processed (data validation failed)
"""
        
        with open(args.stats, 'w') as f:
            f.write(error_stats)
        
        print(error_stats)
        sys.exit(1)
    
    print(f"Processing {len(df)} rows for GO-KFI calculation...")
    
    # GO-KFIを計算
    kfi_df = calculate_go_kfi(df)
    
    if len(kfi_df) == 0:
        print("WARNING: No GO terms found for KFI calculation!", file=sys.stderr)
        print("This may indicate empty GO annotation data.", file=sys.stderr)
    
    # 結果を保存
    try:
        kfi_df.to_csv(args.output, index=False)
        print(f"GO-KFI table saved to: {args.output}")
        print(f"Generated {len(kfi_df)} GO term entries")
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
