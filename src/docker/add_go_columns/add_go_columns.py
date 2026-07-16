#!/usr/bin/env python3
"""
GO列追加スクリプト
UniProt IDからGO annotationを取得してTSVファイルに追加
"""

import sys
import csv
import os
import requests
import time
import itertools
import re
import tempfile

# CSVフィールドサイズ制限を増加
csv.field_size_limit(1 << 30)

BATCH_SIZE = 100  # バッチサイズ
APIURL = "https://rest.uniprot.org/uniprotkb/stream"
API_MAX_ATTEMPTS = 3
GO_ID_PATTERN = re.compile(r'GO:\d{7}')


class UniProtAPIError(RuntimeError):
    """Raised when a UniProt batch cannot be retrieved reliably."""


def _network_connection_guidance(error):
    """Return actionable guidance for connection and timeout failures."""
    network_error_types = tuple(
        error_type
        for error_type in (
            getattr(requests.exceptions, 'ConnectionError', None),
            getattr(requests.exceptions, 'Timeout', None),
        )
        if isinstance(error_type, type)
    )
    if network_error_types and isinstance(error, network_error_types):
        return (
            " Network access is required for UniProt GO annotation. "
            "Ensure this computer and Docker can reach "
            "https://rest.uniprot.org, connect to the internet if needed, "
            "then rerun the tWhEEL go_analysis step."
        )
    return ''


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
    
    response = None
    last_error = None
    for attempt in range(1, API_MAX_ATTEMPTS + 1):
        try:
            response = requests.get(
                APIURL, params=params, headers=headers, timeout=30
            )
            response.raise_for_status()
            break
        except requests.exceptions.RequestException as error:
            last_error = error
            if attempt == API_MAX_ATTEMPTS:
                raise UniProtAPIError(
                    "UniProt API request failed after "
                    f"{API_MAX_ATTEMPTS} attempts for "
                    f"{len(uniprot_ids)} accession(s): {error}."
                    f"{_network_connection_guidance(error)}"
                ) from error

            delay_seconds = 2 ** (attempt - 1)
            print(
                "  UniProt API request failed; retrying in "
                f"{delay_seconds} second(s): {error}",
                file=sys.stderr,
            )
            time.sleep(delay_seconds)

    if response is None:
        raise UniProtAPIError(
            "UniProt API request failed without a response: "
            f"{last_error}"
        )

    result = {}
    reader = csv.reader(response.text.splitlines(), delimiter='\t')
    response_header = next(reader, None)
    if (
        response_header is None
        or len(response_header) < 4
        or response_header[0].strip().lower() not in {'entry', 'accession'}
    ):
        raise UniProtAPIError(
            "UniProt API returned an invalid TSV response for "
            f"{len(uniprot_ids)} accession(s)"
        )

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

    # A successful response may legitimately omit unknown or unannotated IDs.
    empty_annotation = classify_go_terms({})
    for accession in uniprot_ids:
        result.setdefault(accession, empty_annotation.copy())

    print(f"  Retrieved GO annotations for {len(result)}/{len(uniprot_ids)} UniProt IDs")
    return result


def extract_uniprot_id(subject_id):
    """subject_idからUniProt IDを抽出"""
    # UniProt IDの形式: sp|P12345|GENE_HUMAN
    if '|' in subject_id:
        parts = subject_id.split('|')
        if len(parts) >= 2:
            return parts[1]
    return subject_id


def add_go_columns(input_file, output_file):
    """Add GO columns and publish the output only after every batch succeeds."""
    output_path = os.path.abspath(output_file)
    output_dir = os.path.dirname(output_path)
    output_basename = os.path.basename(output_path)
    temporary_output = None

    try:
        with tempfile.NamedTemporaryFile(
            mode='w',
            newline='',
            dir=output_dir,
            prefix=f'.{output_basename}.',
            suffix='.tmp',
            delete=False,
        ) as fout:
            temporary_output = fout.name

            with open(input_file, 'r', newline='') as fin:
                reader = csv.reader(fin, delimiter='\t')
                writer = csv.writer(fout, delimiter='\t')

                # Process the header before requesting any remote data.
                try:
                    header = next(reader)
                except StopIteration as error:
                    raise ValueError("Input file is empty") from error
                new_header = header + ['GO_BP', 'GO_MF', 'GO_CC', 'GO_ALL']
                writer.writerow(new_header)

                go_cache = {}
                processed_rows = 0
                print(
                    "UniProt REST API access is required to retrieve GO "
                    f"annotations: {APIURL}",
                    flush=True,
                )

                while True:
                    chunk = list(itertools.islice(reader, BATCH_SIZE))
                    if not chunk:
                        break

                    uniprot_ids = []
                    for row in chunk:
                        if len(row) > 1:
                            uniprot_ids.append(extract_uniprot_id(row[1]))
                        else:
                            uniprot_ids.append('')

                    new_ids = [
                        uid for uid in uniprot_ids
                        if uid and uid not in go_cache
                    ]
                    if new_ids:
                        # A failed batch raises and prevents publication of a
                        # partially annotated table.
                        go_cache.update(fetch_go_batch(new_ids))
                        time.sleep(1)

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
                            new_row = row + ['', '', '', '']

                        writer.writerow(new_row)
                        processed_rows += 1

                    if processed_rows % 100 == 0:
                        print(f"Processed {processed_rows} rows...")

        os.chmod(temporary_output, 0o644)
        os.replace(temporary_output, output_path)
        temporary_output = None
        return processed_rows
    finally:
        if temporary_output and os.path.exists(temporary_output):
            os.unlink(temporary_output)


def main():
    if len(sys.argv) != 3:
        print("Usage: python add_go_columns.py input.tsv output.tsv", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    print(f"Reading input file: {input_file}", flush=True)

    try:
        processed_rows = add_go_columns(input_file, output_file)
        print(
            f"Successfully processed {processed_rows} rows and wrote to "
            f"{output_file}"
        )
    except Exception as e:
        print(f"Error processing files: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
