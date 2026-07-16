#!/usr/bin/env python3
"""Normalize gene annotations and guide targets into versioned TSV artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
import sqlite3
import tempfile
import uuid
from collections import OrderedDict, defaultdict
from pathlib import Path


SCHEMA_VERSION = 1
NORMALIZED_GENES_FILENAME = 'normalized_genes.tsv'
NORMALIZED_GUIDES_FILENAME = 'normalized_guides.tsv'
MANIFEST_FILENAME = 'normalized_gene_guide_manifest.json'
CURRENT_POINTER_FILENAME = 'normalized_gene_guide_current.json'
LEGACY_FINAL_FILENAME = 'annotation_with_guides_final.tsv'
GENERATIONS_DIRECTORY = '.normalized_gene_guide_generations'

GENE_COLUMNS = (
    'gene_ordinal', 'unique_id', 'classification', 'subject_id', 'gene_name',
    'description', 'annotation_source_database', 'sgd_id', 'ensembl_id',
    'representative_query_id', 'query_ids', 'annotation_row_count', 'identity',
    'length', 'evalue', 'bit_score', 'guide_count', 'has_guides',
)

GUIDE_COLUMNS = (
    'guide_ordinal', 'unique_id', 'reference', 'start', 'end',
    'guide_sequence', 'mm0', 'mm1_2', 'off_target_score', 'strand',
    'raw_guide_target', 'parse_status',
)

_REQUIRED_ANNOTATION_COLUMNS = {
    'query_id', 'subject_id', 'classification', 'unique_id',
}
_IDENTITY_FIELDS = (
    'classification', 'subject_id', 'annotation_source_database', 'sgd_id',
    'ensembl_id', 'gene_name', 'description',
)
_ARTIFACTS = {
    'normalized_genes': NORMALIZED_GENES_FILENAME,
    'normalized_guides': NORMALIZED_GUIDES_FILENAME,
    'annotation_with_guides_final': LEGACY_FINAL_FILENAME,
}


class GeneGuideNormalizationError(ValueError):
    """Raised when a source or normalized artifact violates the schema."""


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_metadata(path, row_count):
    path = Path(path)
    return {
        'filename': path.name,
        'size': path.stat().st_size,
        'sha256': _sha256(path),
        'row_count': int(row_count),
    }


def _source_stat(path):
    source_path = Path(path)
    if source_path.is_symlink() or not source_path.is_file():
        raise GeneGuideNormalizationError(
            f'Normalization source is not a regular file: {source_path}'
        )
    path = source_path.resolve(strict=True)
    stat = path.stat()
    return {
        'path': str(path),
        'size': stat.st_size,
        'mtime_ns': stat.st_mtime_ns,
    }


def _source_metadata(path):
    metadata = _source_stat(path)
    metadata['sha256'] = _sha256(metadata['path'])
    return metadata


def _fsync_file(path):
    with Path(path).open('rb') as handle:
        os.fsync(handle.fileno())


def _fsync_directory(path):
    descriptor = os.open(str(path), os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_json_atomic(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=str(path.parent), prefix=f'.{path.name}.', suffix='.tmp'
    )
    try:
        with os.fdopen(descriptor, 'w', encoding='utf-8') as handle:
            json.dump(payload, handle, sort_keys=True, indent=2)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
        _fsync_directory(path.parent)
    finally:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass


def _number(value, default):
    try:
        number = float(str(value or '').strip())
    except (TypeError, ValueError):
        return default
    return number if math.isfinite(number) else default


def _representative_key(row):
    return (
        _number(row.get('evalue'), math.inf),
        -_number(row.get('bit_score'), -math.inf),
        -_number(row.get('identity'), -math.inf),
        -_number(row.get('length'), -math.inf),
        str(row.get('query_id') or ''),
        tuple((key, str(row.get(key) or '')) for key in sorted(row)),
    )


def _read_annotation(annotation_path):
    with Path(annotation_path).open('r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fieldnames = list(reader.fieldnames or [])
        if len(fieldnames) != len(set(fieldnames)):
            raise GeneGuideNormalizationError(
                'Annotation table contains duplicate columns'
            )
        missing = sorted(_REQUIRED_ANNOTATION_COLUMNS - set(fieldnames))
        if missing:
            raise GeneGuideNormalizationError(
                'Annotation table is missing required columns: ' + ','.join(missing)
            )
        if 'guide_targets' in fieldnames:
            raise GeneGuideNormalizationError(
                'Annotation provenance input must not contain guide_targets'
            )
        rows = []
        for row_number, source_row in enumerate(reader, 2):
            if source_row.get(None):
                raise GeneGuideNormalizationError(
                    f'Annotation row {row_number} contains extra fields'
                )
            row = {name: (source_row.get(name) or '').strip() for name in fieldnames}
            if 'length' not in row and 'alignment_length' in row:
                row['length'] = row['alignment_length']
            unique_id = row.get('unique_id', '')
            classification = row.get('classification', '')
            if not unique_id or not row.get('subject_id') or not row.get('query_id'):
                raise GeneGuideNormalizationError(
                    f'Annotation row {row_number} has an empty required identifier'
                )
            if classification not in {'Multi', 'Single'}:
                raise GeneGuideNormalizationError(
                    f'Annotation row {row_number} has unsupported classification: '
                    f'{classification!r}'
                )
            if not unique_id.startswith(classification + '_'):
                raise GeneGuideNormalizationError(
                    f'Annotation row {row_number} has inconsistent unique_id: '
                    f'{unique_id!r}'
                )
            rows.append(row)
    if not rows:
        raise GeneGuideNormalizationError(
            'Annotation table contains no data rows'
        )
    return fieldnames, rows


def _canonical_genes(annotation_rows):
    groups = OrderedDict()
    subject_to_unique_id = {}
    for row in annotation_rows:
        subject_id = row['subject_id']
        previous_unique_id = subject_to_unique_id.get(subject_id)
        if (
            previous_unique_id is not None
            and previous_unique_id != row['unique_id']
        ):
            raise GeneGuideNormalizationError(
                f'Subject {subject_id!r} maps to conflicting unique IDs: '
                f'{previous_unique_id!r} and {row["unique_id"]!r}'
            )
        subject_to_unique_id[subject_id] = row['unique_id']
        groups.setdefault(row['unique_id'], []).append(row)

    genes = []
    for ordinal, (unique_id, rows) in enumerate(groups.items()):
        stable = {}
        for field in _IDENTITY_FIELDS:
            values = sorted({
                str(row.get(field) or '').strip() for row in rows
                if str(row.get(field) or '').strip()
            })
            if len(values) > 1:
                raise GeneGuideNormalizationError(
                    f'Conflicting {field} values for {unique_id}: {values}'
                )
            stable[field] = values[0] if values else ''

        representative = min(rows, key=_representative_key)
        query_ids = sorted({row['query_id'] for row in rows})
        if len(query_ids) != len(rows):
            raise GeneGuideNormalizationError(
                f'Duplicate query rows for {unique_id}'
            )
        if stable['classification'] == 'Single' and len(query_ids) != 1:
            raise GeneGuideNormalizationError(
                f'Single group {unique_id} must contain exactly one query'
            )
        if stable['classification'] == 'Multi' and len(query_ids) < 2:
            raise GeneGuideNormalizationError(
                f'Multi group {unique_id} must contain at least two queries'
            )
        subject_id = stable['subject_id']
        if subject_id.startswith('SGDID:S') and stable['sgd_id']:
            if stable['sgd_id'] != subject_id[len('SGDID:'):]:
                raise GeneGuideNormalizationError(
                    f'Conflicting subject_id and sgd_id for {unique_id}'
                )
        genes.append({
            'gene_ordinal': ordinal,
            'unique_id': unique_id,
            'classification': stable['classification'],
            'subject_id': stable['subject_id'],
            'gene_name': stable['gene_name'],
            'description': stable['description'],
            'annotation_source_database': stable['annotation_source_database'],
            'sgd_id': stable['sgd_id'],
            'ensembl_id': stable['ensembl_id'],
            'representative_query_id': representative['query_id'],
            'query_ids': query_ids,
            'annotation_row_count': len(rows),
            'identity': representative.get('identity', ''),
            'length': representative.get('length', ''),
            'evalue': representative.get('evalue', ''),
            'bit_score': representative.get('bit_score', ''),
            'guide_count': 0,
            'has_guides': False,
        })
    return genes


def _optional_int(value, field, raw_token):
    if value == '':
        return None
    try:
        return int(value)
    except ValueError as error:
        raise GeneGuideNormalizationError(
            f'Invalid {field} in guide target {raw_token!r}'
        ) from error


def _optional_float(value, field, raw_token):
    if value == '':
        return None
    try:
        number = float(value)
    except ValueError as error:
        raise GeneGuideNormalizationError(
            f'Invalid {field} in guide target {raw_token!r}'
        ) from error
    if not math.isfinite(number):
        raise GeneGuideNormalizationError(
            f'Invalid {field} in guide target {raw_token!r}'
        )
    return number


def _parse_guide_token(raw_token, expected_classification, known_gene_ids):
    raw_token = raw_token.strip()
    if not raw_token:
        raise GeneGuideNormalizationError('Guide target token is empty')
    parts = [part.strip() for part in raw_token.split('/')]
    if len(parts) == 1:
        return {
            'unique_id': '', 'reference': '', 'start': None, 'end': None,
            'guide_sequence': '', 'mm0': None, 'mm1_2': None,
            'off_target_score': None, 'strand': '',
            'raw_guide_target': raw_token, 'parse_status': 'opaque_legacy',
        }
    if len(parts) not in {6, 7, 8}:
        raise GeneGuideNormalizationError(
            f'Guide target must have 6, 7, or 8 parts: {raw_token!r}'
        )

    reference_id = parts[0]
    if not reference_id or not parts[3]:
        raise GeneGuideNormalizationError(
            f'Guide reference and sequence must not be empty: {raw_token!r}'
        )
    unique_id = reference_id
    suffix = '_aligned_consensus'
    if expected_classification == 'Multi' and unique_id.endswith(suffix):
        unique_id = unique_id[:-len(suffix)]
    if unique_id not in known_gene_ids:
        raise GeneGuideNormalizationError(
            f'Guide target references unknown gene {unique_id!r}'
        )
    if known_gene_ids[unique_id] != expected_classification:
        raise GeneGuideNormalizationError(
            f'Guide target classification does not match {unique_id}'
        )

    start = _optional_int(parts[1], 'start', raw_token)
    end = _optional_int(parts[2], 'end', raw_token)
    mm0 = _optional_int(parts[4], 'mm0', raw_token)
    mm1_2 = _optional_int(parts[5], 'mm1_2', raw_token)
    if start is None or end is None or start < 0 or end <= start:
        raise GeneGuideNormalizationError(
            f'Invalid guide coordinates in {raw_token!r}'
        )
    if mm0 is None or mm1_2 is None or mm0 < 0 or mm1_2 < 0:
        raise GeneGuideNormalizationError(
            f'Invalid guide mismatch counts in {raw_token!r}'
        )
    strand = parts[7] if len(parts) == 8 else ''
    score_text = parts[6] if len(parts) >= 7 else ''
    if len(parts) == 7 and parts[6] in {'+', '-'}:
        strand = parts[6]
        score_text = ''
    if strand not in {'', '+', '-'}:
        raise GeneGuideNormalizationError(
            f'Invalid strand in guide target {raw_token!r}'
        )
    return {
        'unique_id': unique_id,
        'reference': reference_id,
        'start': start,
        'end': end,
        'guide_sequence': parts[3].upper(),
        'mm0': mm0,
        'mm1_2': mm1_2,
        'off_target_score': (
            _optional_float(score_text, 'off_target_score', raw_token)
            if score_text else None
        ),
        'strand': strand,
        'raw_guide_target': raw_token,
        'parse_status': f'parsed_{len(parts)}_part',
    }


def _guide_tokens(path):
    if path is None:
        return
    with Path(path).open('r', encoding='utf-8', newline='') as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.rstrip('\r\n')
            if not stripped.strip():
                continue
            fields = [field.strip() for field in stripped.split('\t')]
            if len(fields) == 1:
                yield fields[0], line_number
            elif len(fields) >= 6:
                yield '/'.join(fields[:6]), line_number
            else:
                raise GeneGuideNormalizationError(
                    f'Malformed guide BED row at {path}:{line_number}'
                )


def _guide_identity(guide):
    """Return the typed guide identity used for de-duplication."""
    identity = (
        guide['unique_id'],
        guide['reference'],
        guide['start'],
        guide['end'],
        guide['guide_sequence'],
        guide['mm0'],
        guide['mm1_2'],
        guide['off_target_score'],
        guide['strand'],
    )
    if guide['parse_status'] == 'opaque_legacy':
        return identity + (guide['raw_guide_target'],)
    return identity


def _create_guide_index(database_path, multi_bed_path, single_bed_path, genes):
    """Build a disk-backed typed guide index with bounded Python memory."""
    known_gene_ids = {
        gene['unique_id']: gene['classification'] for gene in genes
    }
    opaque_gene_ids = [
        gene['unique_id'] for gene in genes
        if gene['classification'] == 'Single'
    ]
    connection = sqlite3.connect(str(database_path))
    try:
        connection.execute('PRAGMA journal_mode=OFF')
        connection.execute('PRAGMA synchronous=OFF')
        connection.execute('PRAGMA temp_store=FILE')
        connection.execute(
            """
            CREATE TABLE guides (
                identity_key TEXT PRIMARY KEY,
                first_ordinal INTEGER NOT NULL,
                unique_id TEXT NOT NULL,
                reference TEXT NOT NULL,
                start INTEGER,
                end INTEGER,
                guide_sequence TEXT NOT NULL,
                mm0 INTEGER,
                mm1_2 INTEGER,
                off_target_score TEXT NOT NULL,
                strand TEXT NOT NULL,
                raw_guide_target TEXT NOT NULL,
                parse_status TEXT NOT NULL
            ) WITHOUT ROWID
            """
        )
        source_ordinal = 0
        for classification, path in (
            ('Multi', multi_bed_path), ('Single', single_bed_path)
        ):
            for raw_token, line_number in _guide_tokens(path) or ():
                guide = _parse_guide_token(
                    raw_token, classification, known_gene_ids
                )
                if guide['parse_status'] == 'opaque_legacy':
                    if len(opaque_gene_ids) != 1:
                        raise GeneGuideNormalizationError(
                            'Opaque legacy guide targets require exactly one '
                            f'Single gene (at {path}:{line_number})'
                        )
                    guide['unique_id'] = opaque_gene_ids[0]
                identity_key = hashlib.sha256(
                    json.dumps(
                        _guide_identity(guide),
                        ensure_ascii=False,
                        separators=(',', ':'),
                    ).encode('utf-8')
                ).hexdigest()
                score = guide['off_target_score']
                score_text = '' if score is None else format(score, '.17g')
                connection.execute(
                    """
                    INSERT OR IGNORE INTO guides (
                        identity_key, first_ordinal, unique_id, reference,
                        start, end, guide_sequence, mm0, mm1_2,
                        off_target_score, strand, raw_guide_target,
                        parse_status
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        identity_key,
                        source_ordinal,
                        guide['unique_id'],
                        guide['reference'],
                        guide['start'],
                        guide['end'],
                        guide['guide_sequence'],
                        guide['mm0'],
                        guide['mm1_2'],
                        score_text,
                        guide['strand'],
                        guide['raw_guide_target'],
                        guide['parse_status'],
                    ),
                )
                source_ordinal += 1
        connection.execute(
            'CREATE UNIQUE INDEX guides_by_ordinal '
            'ON guides(first_ordinal)'
        )
        connection.execute(
            'CREATE INDEX guides_by_gene_ordinal '
            'ON guides(unique_id, first_ordinal)'
        )
        connection.commit()
        return connection
    except Exception:
        connection.close()
        raise


def _guide_rows_from_index(connection, unique_id=None):
    columns = (
        'unique_id', 'reference', 'start', 'end', 'guide_sequence', 'mm0',
        'mm1_2', 'off_target_score', 'strand', 'raw_guide_target',
        'parse_status',
    )
    query = (
        'SELECT ' + ','.join(columns) + ' FROM guides'
    )
    parameters = ()
    if unique_id is not None:
        query += ' WHERE unique_id = ?'
        parameters = (unique_id,)
    query += ' ORDER BY first_ordinal'
    for ordinal, values in enumerate(connection.execute(query, parameters)):
        row = dict(zip(columns, values))
        row['guide_ordinal'] = ordinal
        yield row


def _guide_counts_from_index(connection):
    return defaultdict(int, {
        unique_id: int(count)
        for unique_id, count in connection.execute(
            'SELECT unique_id, COUNT(*) FROM guides GROUP BY unique_id'
        )
    })


def _write_guides_from_index(path, connection):
    count = 0
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(
            handle, fieldnames=GUIDE_COLUMNS, delimiter='\t',
            lineterminator='\n',
        )
        writer.writeheader()
        for row in _guide_rows_from_index(connection):
            writer.writerow(row)
            count += 1
        handle.flush()
        os.fsync(handle.fileno())
    return count


def _write_tsv(path, columns, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        for row in rows:
            serialized = dict(row)
            if 'query_ids' in serialized:
                serialized['query_ids'] = json.dumps(
                    serialized['query_ids'], ensure_ascii=False, separators=(',', ':')
                )
            for key, value in list(serialized.items()):
                if value is None:
                    serialized[key] = ''
                elif isinstance(value, bool):
                    serialized[key] = 'true' if value else 'false'
            writer.writerow(serialized)
        handle.flush()
        os.fsync(handle.fileno())


def _write_legacy_final(path, fieldnames, annotation_rows, guide_connection):
    """Write a compact one-row-per-gene compatibility table."""
    annotation_groups = OrderedDict()
    for row in annotation_rows:
        annotation_groups.setdefault(row['unique_id'], []).append(row)
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[*fieldnames, 'guide_targets'],
            delimiter='\t',
            lineterminator='\n',
        )
        writer.writeheader()
        for unique_id, rows in annotation_groups.items():
            representative = min(rows, key=_representative_key)
            output = {
                name: representative.get(name, '') for name in fieldnames
            }
            output['guide_targets'] = ','.join(
                guide['raw_guide_target']
                for guide in _guide_rows_from_index(
                    guide_connection, unique_id=unique_id
                )
            )
            writer.writerow(output)
        handle.flush()
        os.fsync(handle.fileno())


def publish_normalized_gene_guide_generation(
    annotation_path,
    output_dir,
    multi_bed_path=None,
    single_bed_path=None,
):
    """Build, validate, and atomically select one normalized generation."""
    output_dir = Path(output_dir)
    generations_dir = output_dir / GENERATIONS_DIRECTORY
    generations_dir.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(tempfile.mkdtemp(prefix='.staging-', dir=generations_dir))
    try:
        source_paths = OrderedDict([
            ('annotation_classified', annotation_path),
            ('multi_guides', multi_bed_path),
            ('single_guides', single_bed_path),
        ])
        source_paths = OrderedDict(
            (key, path) for key, path in source_paths.items()
            if path is not None
        )
        source_stats_before = {
            key: _source_stat(path) for key, path in source_paths.items()
        }
        fieldnames, annotation_rows = _read_annotation(annotation_path)
        genes = _canonical_genes(annotation_rows)
        genes_path = staging_dir / NORMALIZED_GENES_FILENAME
        guides_path = staging_dir / NORMALIZED_GUIDES_FILENAME
        legacy_path = staging_dir / LEGACY_FINAL_FILENAME
        manifest_path = staging_dir / MANIFEST_FILENAME
        guide_index_path = staging_dir / '.guide_index.sqlite3'
        guide_connection = _create_guide_index(
            guide_index_path,
            multi_bed_path,
            single_bed_path,
            genes,
        )
        try:
            guide_counts = _guide_counts_from_index(guide_connection)
            for gene in genes:
                gene['guide_count'] = guide_counts[gene['unique_id']]
                gene['has_guides'] = bool(gene['guide_count'])

            _write_tsv(genes_path, GENE_COLUMNS, genes)
            guide_count = _write_guides_from_index(
                guides_path, guide_connection
            )
            _write_legacy_final(
                legacy_path,
                fieldnames,
                annotation_rows,
                guide_connection,
            )
        finally:
            guide_connection.close()
        guide_index_path.unlink(missing_ok=True)
        source_metadata = {
            key: _source_metadata(path)
            for key, path in source_paths.items()
        }
        for key, previous in source_stats_before.items():
            current = source_metadata[key]
            if any(
                current[field] != previous[field]
                for field in ('path', 'size', 'mtime_ns')
            ):
                raise GeneGuideNormalizationError(
                    f'Normalization source changed during publication: {key}'
                )
        manifest = {
            'schema_version': SCHEMA_VERSION,
            'sources': source_metadata,
            'artifacts': {
                'normalized_genes': _artifact_metadata(genes_path, len(genes)),
                'normalized_guides': _artifact_metadata(
                    guides_path, guide_count
                ),
                'annotation_with_guides_final': _artifact_metadata(
                    legacy_path, len(genes)
                ),
            },
            'gene_count': len(genes),
            'guide_count': guide_count,
            'annotation_row_count': len(annotation_rows),
            'legacy_row_count': len(genes),
        }
        _write_json_atomic(manifest_path, manifest)
        validate_normalized_gene_guide_generation(
            genes_path,
            guides_path,
            manifest_path,
            legacy_path=legacy_path,
        )
        generation_id = uuid.uuid4().hex
        generation_dir = generations_dir / generation_id
        os.replace(staging_dir, generation_dir)
        _fsync_directory(generations_dir)
        staging_dir = None
        pointer_path = output_dir / CURRENT_POINTER_FILENAME
        _write_json_atomic(pointer_path, {
            'schema_version': SCHEMA_VERSION,
            'generation': generation_id,
            'manifest': MANIFEST_FILENAME,
            'manifest_size': (
                generation_dir / MANIFEST_FILENAME
            ).stat().st_size,
            'manifest_sha256': _sha256(
                generation_dir / MANIFEST_FILENAME
            ),
        })
        resolved = resolve_current_normalized_gene_guide_paths(output_dir)
        return {
            **resolved,
            'annotation_with_guides_final': str(
                generation_dir / LEGACY_FINAL_FILENAME
            ),
            'normalized_gene_guide_current': str(pointer_path),
            'generation_dir': str(generation_dir),
        }
    finally:
        if staging_dir is not None:
            shutil.rmtree(staging_dir, ignore_errors=True)


def _load_manifest(manifest_path):
    try:
        with Path(manifest_path).open('r', encoding='utf-8') as handle:
            manifest = json.load(handle)
    except (OSError, ValueError, TypeError) as error:
        raise GeneGuideNormalizationError(
            f'Cannot read normalized gene-guide manifest: {manifest_path}'
        ) from error
    if manifest.get('schema_version') != SCHEMA_VERSION:
        raise GeneGuideNormalizationError('Unsupported normalized schema version')
    return manifest


def _verify_artifact(path, metadata, logical_key):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise GeneGuideNormalizationError(
            f'Normalized artifact is missing: {logical_key}'
        )
    if not isinstance(metadata, dict):
        raise GeneGuideNormalizationError(
            f'Normalized artifact metadata is invalid: {logical_key}'
        )
    if path.name != metadata.get('filename'):
        raise GeneGuideNormalizationError(
            f'Normalized artifact filename mismatch: {logical_key}'
        )
    expected_size = metadata.get('size')
    expected_hash = metadata.get('sha256')
    expected_rows = metadata.get('row_count')
    if (
        not isinstance(expected_size, int)
        or expected_size < 0
        or not isinstance(expected_rows, int)
        or expected_rows < 0
        or not isinstance(expected_hash, str)
        or re.fullmatch(r'[0-9a-f]{64}', expected_hash) is None
    ):
        raise GeneGuideNormalizationError(
            f'Normalized artifact metadata is invalid: {logical_key}'
        )
    if path.stat().st_size != expected_size:
        raise GeneGuideNormalizationError(
            f'Normalized artifact size mismatch: {logical_key}'
        )
    if _sha256(path) != expected_hash:
        raise GeneGuideNormalizationError(
            f'Normalized artifact hash mismatch: {logical_key}'
        )


def _verify_source_metadata(metadata, logical_key):
    if not isinstance(metadata, dict):
        raise GeneGuideNormalizationError(
            f'Normalization source metadata is invalid: {logical_key}'
        )
    path_value = metadata.get('path')
    expected_size = metadata.get('size')
    expected_mtime = metadata.get('mtime_ns')
    expected_hash = metadata.get('sha256')
    if (
        not isinstance(path_value, str)
        or not isinstance(expected_size, int)
        or expected_size < 0
        or not isinstance(expected_mtime, int)
        or expected_mtime < 0
        or not isinstance(expected_hash, str)
        or re.fullmatch(r'[0-9a-f]{64}', expected_hash) is None
    ):
        raise GeneGuideNormalizationError(
            f'Normalization source metadata is invalid: {logical_key}'
        )
    current = _source_stat(path_value)
    if (
        current['path'] != path_value
        or current['size'] != expected_size
        or current['mtime_ns'] != expected_mtime
        or _sha256(path_value) != expected_hash
    ):
        raise GeneGuideNormalizationError(
            f'Normalized generation source is stale: {logical_key}'
        )


def resolve_current_normalized_gene_guide_paths(output_dir):
    """Resolve and integrity-check current artifacts without loading TSV rows."""
    output_dir = Path(output_dir)
    pointer_path = output_dir / CURRENT_POINTER_FILENAME
    if not os.path.lexists(str(pointer_path)):
        return {}
    if pointer_path.is_symlink() or not pointer_path.is_file():
        raise GeneGuideNormalizationError(
            f'Invalid normalized generation pointer: {pointer_path}'
        )
    try:
        with pointer_path.open('r', encoding='utf-8') as handle:
            pointer = json.load(handle)
    except (OSError, ValueError, TypeError) as error:
        raise GeneGuideNormalizationError(
            f'Cannot read normalized generation pointer: {pointer_path}'
        ) from error
    if pointer.get('schema_version') != SCHEMA_VERSION:
        raise GeneGuideNormalizationError('Unsupported pointer schema version')
    generation = pointer.get('generation')
    if (
        not isinstance(generation, str)
        or re.fullmatch(r'[0-9a-f]{32}', generation) is None
    ):
        raise GeneGuideNormalizationError('Invalid normalized generation pointer')
    generations_dir = (output_dir / GENERATIONS_DIRECTORY).resolve()
    generation_dir = generations_dir / generation
    try:
        resolved_generation_dir = generation_dir.resolve(strict=True)
        resolved_generation_dir.relative_to(generations_dir)
    except (OSError, ValueError) as error:
        raise GeneGuideNormalizationError(
            'Normalized generation path is invalid'
        ) from error
    if generation_dir.is_symlink() or not resolved_generation_dir.is_dir():
        raise GeneGuideNormalizationError(
            'Normalized generation path is invalid'
        )
    generation_dir = resolved_generation_dir
    if pointer.get('manifest', MANIFEST_FILENAME) != MANIFEST_FILENAME:
        raise GeneGuideNormalizationError('Invalid normalized manifest pointer')
    manifest_path = generation_dir / MANIFEST_FILENAME
    if manifest_path.is_symlink() or not manifest_path.is_file():
        raise GeneGuideNormalizationError('Normalized manifest is missing')
    if manifest_path.stat().st_size != pointer.get('manifest_size'):
        raise GeneGuideNormalizationError(
            'Normalized current manifest size mismatch'
        )
    if _sha256(manifest_path) != pointer.get('manifest_sha256'):
        raise GeneGuideNormalizationError(
            'Normalized current manifest hash mismatch'
        )
    manifest = _load_manifest(manifest_path)
    sources = manifest.get('sources')
    if (
        not isinstance(sources, dict)
        or 'annotation_classified' not in sources
    ):
        raise GeneGuideNormalizationError(
            'Normalized manifest source metadata is missing'
        )
    for logical_key, metadata in sorted(sources.items()):
        _verify_source_metadata(metadata, logical_key)
    artifacts = manifest.get('artifacts')
    if not isinstance(artifacts, dict):
        raise GeneGuideNormalizationError('Manifest artifacts are missing')
    for logical_key in ('normalized_genes', 'normalized_guides'):
        filename = _ARTIFACTS[logical_key]
        metadata = artifacts.get(logical_key)
        if not isinstance(metadata, dict):
            raise GeneGuideNormalizationError(
                f'Manifest artifact metadata is missing: {logical_key}'
            )
        _verify_artifact(generation_dir / filename, metadata, logical_key)
    return {
        'normalized_genes': str(generation_dir / NORMALIZED_GENES_FILENAME),
        'normalized_guides': str(generation_dir / NORMALIZED_GUIDES_FILENAME),
        'normalized_gene_guide_manifest': str(manifest_path),
    }


def _parse_required_int(value, field, row_number):
    try:
        return int(value)
    except (TypeError, ValueError) as error:
        raise GeneGuideNormalizationError(
            f'Invalid {field} at normalized row {row_number}'
        ) from error


def _parse_bool(value, field, row_number):
    normalized = str(value).strip().lower()
    if normalized == 'true':
        return True
    if normalized == 'false':
        return False
    raise GeneGuideNormalizationError(
        f'Invalid {field} at normalized row {row_number}'
    )


def _load_typed_genes(genes_path):
    genes = []
    with Path(genes_path).open('r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if tuple(reader.fieldnames or ()) != GENE_COLUMNS:
            raise GeneGuideNormalizationError('Normalized gene columns do not match')
        for row_number, row in enumerate(reader, 2):
            typed = dict(row)
            typed['gene_ordinal'] = _parse_required_int(
                row['gene_ordinal'], 'gene_ordinal', row_number
            )
            typed['annotation_row_count'] = _parse_required_int(
                row['annotation_row_count'], 'annotation_row_count', row_number
            )
            typed['guide_count'] = _parse_required_int(
                row['guide_count'], 'guide_count', row_number
            )
            typed['has_guides'] = _parse_bool(
                row['has_guides'], 'has_guides', row_number
            )
            if typed['gene_ordinal'] != len(genes):
                raise GeneGuideNormalizationError(
                    f'Non-contiguous gene_ordinal at normalized row '
                    f'{row_number}'
                )
            if typed['annotation_row_count'] <= 0 or typed['guide_count'] < 0:
                raise GeneGuideNormalizationError(
                    f'Invalid normalized gene counts at row {row_number}'
                )
            if typed['has_guides'] != (typed['guide_count'] > 0):
                raise GeneGuideNormalizationError(
                    f'has_guides does not match guide_count at row '
                    f'{row_number}'
                )
            try:
                typed['query_ids'] = json.loads(row['query_ids'])
            except (ValueError, TypeError) as error:
                raise GeneGuideNormalizationError(
                    f'Invalid query_ids at normalized row {row_number}'
                ) from error
            if (
                not isinstance(typed['query_ids'], list)
                or not typed['query_ids']
                or not all(
                    isinstance(value, str) and value
                    for value in typed['query_ids']
                )
                or typed['query_ids'] != sorted(set(typed['query_ids']))
                or typed['annotation_row_count'] != len(typed['query_ids'])
            ):
                raise GeneGuideNormalizationError(
                    f'Invalid query_ids at normalized row {row_number}'
                )
            genes.append(typed)
    return genes


def _iter_typed_guides(guides_path):
    with Path(guides_path).open('r', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if tuple(reader.fieldnames or ()) != GUIDE_COLUMNS:
            raise GeneGuideNormalizationError('Normalized guide columns do not match')
        for row_number, row in enumerate(reader, 2):
            typed = dict(row)
            typed['guide_ordinal'] = _parse_required_int(
                row['guide_ordinal'], 'guide_ordinal', row_number
            )
            if typed['guide_ordinal'] != row_number - 2:
                raise GeneGuideNormalizationError(
                    f'Non-contiguous guide_ordinal at normalized row '
                    f'{row_number}'
                )
            for field in ('start', 'end', 'mm0', 'mm1_2'):
                typed[field] = (
                    None if row[field] == '' else _parse_required_int(
                        row[field], field, row_number
                    )
                )
            typed['off_target_score'] = (
                None if row['off_target_score'] == '' else
                _optional_float(
                    row['off_target_score'], 'off_target_score',
                    row['raw_guide_target'],
                )
            )
            if typed['parse_status'] not in {
                'parsed_6_part', 'parsed_7_part', 'parsed_8_part',
                'opaque_legacy',
            }:
                raise GeneGuideNormalizationError(
                    f'Invalid parse_status at normalized row {row_number}'
                )
            if typed['parse_status'] != 'opaque_legacy':
                if (
                    not typed['reference']
                    or not typed['guide_sequence']
                    or any(
                        typed[field] is None or typed[field] < 0
                        for field in ('start', 'end', 'mm0', 'mm1_2')
                    )
                    or typed['end'] <= typed['start']
                    or typed['strand'] not in {'', '+', '-'}
                ):
                    raise GeneGuideNormalizationError(
                        f'Invalid typed guide at normalized row {row_number}'
                    )
            yield typed


def load_normalized_gene_data(
    genes_path, guides_path, manifest_path, verify=True
):
    """Load typed genes after validating both normalized artifact hashes."""
    manifest = _load_manifest(manifest_path)
    artifact_paths = {
        'normalized_genes': Path(genes_path),
        'normalized_guides': Path(guides_path),
    }
    for key, path in artifact_paths.items():
        metadata = (manifest.get('artifacts') or {}).get(key)
        if not isinstance(metadata, dict):
            raise GeneGuideNormalizationError(
                f'Manifest artifact metadata is missing: {key}'
            )
        if verify:
            _verify_artifact(path, metadata, key)
        elif not path.is_file():
            raise GeneGuideNormalizationError(
                f'Normalized artifact is missing: {key}'
            )
    genes = _load_typed_genes(genes_path)
    expected_gene_count = manifest.get('gene_count')
    expected_rows = manifest['artifacts']['normalized_genes']['row_count']
    if len(genes) != expected_gene_count or len(genes) != expected_rows:
        raise GeneGuideNormalizationError(
            'Normalized gene row counts do not match manifest'
        )
    gene_ids = [gene['unique_id'] for gene in genes]
    if (
        any(
            not unique_id or gene['classification'] not in {'Multi', 'Single'}
            or not gene['subject_id']
            for unique_id, gene in zip(gene_ids, genes)
        )
        or len(gene_ids) != len(set(gene_ids))
    ):
        raise GeneGuideNormalizationError(
            'Normalized genes contain invalid or duplicate identities'
        )
    return {
        'genes': genes,
        'manifest': manifest,
        'paths': {
            'normalized_genes': str(genes_path),
            'normalized_guides': str(guides_path),
            'normalized_gene_guide_manifest': str(manifest_path),
        },
    }


def iter_normalized_guides(guides_path, manifest_path, verify=True):
    """Yield typed guides without materializing the complete guide table."""
    manifest = _load_manifest(manifest_path)
    metadata = (manifest.get('artifacts') or {}).get('normalized_guides')
    if not isinstance(metadata, dict):
        raise GeneGuideNormalizationError(
            'Manifest artifact metadata is missing: normalized_guides'
        )
    if verify:
        _verify_artifact(guides_path, metadata, 'normalized_guides')
    observed_rows = 0
    for guide in _iter_typed_guides(guides_path):
        observed_rows += 1
        yield guide
    if (
        observed_rows != metadata.get('row_count')
        or observed_rows != manifest.get('guide_count')
    ):
        raise GeneGuideNormalizationError(
            'Normalized guide row counts do not match manifest'
        )


def validate_normalized_gene_guide_generation(
    genes_path, guides_path, manifest_path, legacy_path=None
):
    """Validate a complete generation with bounded guide-table memory."""
    dataset = load_normalized_gene_data(
        genes_path, guides_path, manifest_path
    )
    genes_by_id = {
        gene['unique_id']: gene for gene in dataset['genes']
    }
    observed_counts = {unique_id: 0 for unique_id in genes_by_id}
    for guide in iter_normalized_guides(
        guides_path, manifest_path, verify=False
    ):
        unique_id = guide['unique_id']
        if unique_id not in genes_by_id:
            raise GeneGuideNormalizationError(
                f'Guide references unknown gene: {unique_id}'
            )
        observed_counts[unique_id] += 1
    for unique_id, gene in genes_by_id.items():
        if observed_counts[unique_id] != gene['guide_count']:
            raise GeneGuideNormalizationError(
                f'Guide count does not reconcile for {unique_id}'
            )

    if legacy_path is not None:
        manifest = dataset['manifest']
        metadata = (manifest.get('artifacts') or {}).get(
            'annotation_with_guides_final'
        )
        _verify_artifact(
            legacy_path, metadata, 'annotation_with_guides_final'
        )
        with Path(legacy_path).open(
            'r', encoding='utf-8', newline=''
        ) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            fieldnames = set(reader.fieldnames or [])
            if not {'unique_id', 'guide_targets'}.issubset(fieldnames):
                raise GeneGuideNormalizationError(
                    'Legacy compatibility table columns do not match'
                )
            legacy_rows = sum(1 for _ in reader)
        if legacy_rows != metadata.get('row_count'):
            raise GeneGuideNormalizationError(
                'Legacy compatibility row count does not match manifest'
            )
    return dataset


def load_normalized_gene_guide_data(
    genes_path,
    guides_path,
    manifest_path,
    *,
    legacy_path=None,
):
    """Load a normalized generation with typed values and relational checks."""
    dataset = load_normalized_gene_data(
        genes_path, guides_path, manifest_path
    )
    manifest = dataset['manifest']
    if legacy_path is not None:
        metadata = (manifest.get('artifacts') or {}).get(
            'annotation_with_guides_final'
        )
        if not isinstance(metadata, dict):
            raise GeneGuideNormalizationError(
                'Manifest artifact metadata is missing: '
                'annotation_with_guides_final'
            )
        _verify_artifact(
            legacy_path, metadata, 'annotation_with_guides_final'
        )

    genes = dataset['genes']
    guides = list(iter_normalized_guides(
        guides_path, manifest_path, verify=False
    ))

    expected_gene_count = manifest.get('gene_count')
    expected_guide_count = manifest.get('guide_count')
    artifacts = manifest.get('artifacts') or {}
    if (
        len(genes) != expected_gene_count
        or len(guides) != expected_guide_count
        or len(genes) != artifacts['normalized_genes'].get('row_count')
        or len(guides) != artifacts['normalized_guides'].get('row_count')
    ):
        raise GeneGuideNormalizationError('Normalized row counts do not match manifest')
    gene_ids = [gene['unique_id'] for gene in genes]
    gene_id_set = set(gene_ids)
    guides_by_unique_id = defaultdict(list)
    identities = set()
    for guide in guides:
        if guide['unique_id'] not in gene_id_set:
            raise GeneGuideNormalizationError(
                f"Guide references unknown gene: {guide['unique_id']}"
            )
        identity = _guide_identity(guide)
        if identity in identities:
            raise GeneGuideNormalizationError('Normalized guides contain duplicates')
        identities.add(identity)
        guides_by_unique_id[guide['unique_id']].append(guide)
    for gene in genes:
        count = len(guides_by_unique_id.get(gene['unique_id'], ()))
        if gene['guide_count'] != count or gene['has_guides'] != bool(count):
            raise GeneGuideNormalizationError(
                f"Guide count does not reconcile for {gene['unique_id']}"
            )
    return {
        'genes': genes,
        'guides': guides,
        'guides_by_unique_id': dict(guides_by_unique_id),
        'manifest': manifest,
    }


def load_current_normalized_gene_guide_data(output_dir):
    """Resolve and fully load the current normalized generation."""
    paths = resolve_current_normalized_gene_guide_paths(output_dir)
    if not paths:
        return {}
    dataset = load_normalized_gene_guide_data(
        paths['normalized_genes'],
        paths['normalized_guides'],
        paths['normalized_gene_guide_manifest'],
    )
    dataset['paths'] = paths
    return dataset
