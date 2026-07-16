#!/usr/bin/env python3
"""
Pathway Report Generation Module
HTML report generation and database integration of PW-KFI analysis results
"""

import os
import json
import logging
import sqlite3
import csv
import ast
import shutil
import html
import re
import sys
import math
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import datetime
from pathlib import Path


try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2 ** 31 - 1)


class PathwayReportGenerator:
    """Pathway report generation processing class"""
    
    def __init__(self, config):
        """
        Initialization
        
        Args:
            config (dict): Configuration Dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        self.output_dir = self.config['output']['base_dir']
        self.web_data_dir = os.path.join(self.output_dir, 'web_data')
        self.pathway_report_config = self.config.get('pathway_report', {})
        
        # 出力ディレクトリを作成
        Path(self.web_data_dir).mkdir(parents=True, exist_ok=True)
    
    def load_pathway_analysis_results(self):
        """Load and validate canonical-gene pathway analysis artifacts."""
        kfi_table_path = os.path.join(self.intermediate_dir, 'pathway_kfi_table.csv')
        if not os.path.exists(kfi_table_path):
            raise FileNotFoundError(f"PW-KFI table not found.: {kfi_table_path}")

        required_pathway_columns = {
            'pathway_id', 'pathway_name', 'total_gene_nodes',
            'evaluable_gene_node_count',
            'unique_xref_count', 'total_genes', 'gene_count_basis',
            'canonicalization_complete', 'canonicalization_status',
            'matched_genes', 'unmatched_genes', 'unmatchable_genes',
            'unmatched_gene_count', 'matched_gene_count',
            'genes_with_guides', 'pw_kfi_score', 'coverage_rate',
            'unmatchable_gene_count',
        }
        pathways = []
        expected_mapping_by_pathway = {}
        expected_pathway_names = {}
        with open(kfi_table_path, 'r', encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f)
            missing_columns = required_pathway_columns - set(
                reader.fieldnames or []
            )
            if missing_columns:
                raise ValueError(
                    "Legacy or incomplete PW-KFI table cannot be reported as "
                    "canonical-gene results. Regenerate pathway analysis; "
                    "missing columns: " + ', '.join(sorted(missing_columns))
                )
            for row in reader:
                pathway_id = row.get('pathway_id', '')
                if not pathway_id:
                    raise ValueError('PW-KFI table contains an empty pathway_id')
                if pathway_id in expected_mapping_by_pathway:
                    raise ValueError(
                        f"PW-KFI table contains duplicate pathway_id: {pathway_id}"
                    )
                if row.get('gene_count_basis') != 'exact_canonical_gene_units':
                    raise ValueError(
                        f"Unsupported gene_count_basis for {pathway_id}: "
                        f"{row.get('gene_count_basis')!r}. Regenerate pathway analysis."
                    )
                if (
                    str(row.get('canonicalization_complete', '')).lower()
                    != 'true'
                    or row.get('canonicalization_status') != 'complete'
                ):
                    raise ValueError(
                        f"PW-KFI is not canonically resolved for {pathway_id}"
                    )

                counts = {
                    field: self._nonnegative_int(
                        row.get(field, ''), field, pathway_id
                    )
                    for field in (
                        'total_gene_nodes', 'evaluable_gene_node_count',
                        'unique_xref_count', 'total_genes',
                        'unmatchable_gene_count', 'matched_gene_count',
                        'unmatched_gene_count', 'genes_with_guides',
                    )
                }
                if not (
                    counts['total_genes'] <= counts['unique_xref_count']
                    <= counts['evaluable_gene_node_count']
                    <= counts['total_gene_nodes']
                ):
                    raise ValueError(
                        f"Invalid pathway count ordering for {pathway_id}: {counts}"
                    )
                if (
                    counts['evaluable_gene_node_count']
                    + counts['unmatchable_gene_count']
                    != counts['total_gene_nodes']
                ):
                    raise ValueError(
                        f"Raw node counts do not reconcile for {pathway_id}"
                    )
                if (
                    counts['matched_gene_count']
                    + counts['unmatched_gene_count']
                    != counts['total_genes']
                ):
                    raise ValueError(
                        f"Canonical gene counts do not reconcile for {pathway_id}"
                    )
                if counts['genes_with_guides'] > counts['matched_gene_count']:
                    raise ValueError(
                        f"Guide-positive genes exceed matched genes for {pathway_id}"
                    )

                serialized_gene_fields = {}
                for field in (
                    'matched_genes', 'unmatched_genes', 'unmatchable_genes'
                ):
                    try:
                        value = self._parse_matched_genes(row.get(field, ''))
                    except (
                        TypeError, ValueError, SyntaxError,
                        json.JSONDecodeError,
                    ) as exc:
                        raise ValueError(
                            f"Cannot deserialize {field} for {pathway_id}"
                        ) from exc
                    if not isinstance(value, list):
                        raise ValueError(
                            f"{field} must be a list for {pathway_id}"
                        )
                    serialized_gene_fields[field] = value

                matched_genes = serialized_gene_fields['matched_genes']
                unmatched_genes = serialized_gene_fields['unmatched_genes']
                unmatchable_genes = serialized_gene_fields['unmatchable_genes']
                if len(matched_genes) != counts['matched_gene_count']:
                    raise ValueError(
                        f"matched_genes does not reconcile for {pathway_id}"
                    )
                if len(unmatched_genes) != counts['unmatched_gene_count']:
                    raise ValueError(
                        f"unmatched_genes does not reconcile for {pathway_id}"
                    )
                if len(unmatchable_genes) != counts['unmatchable_gene_count']:
                    raise ValueError(
                        f"unmatchable_genes does not reconcile for {pathway_id}"
                    )
                if any(not isinstance(gene, dict) for gene in unmatchable_genes):
                    raise ValueError(
                        f"unmatchable_genes contains an invalid record for "
                        f"{pathway_id}"
                    )

                canonical_keys = []
                source_xref_keys = []
                evaluable_node_count = 0
                guide_positive_count = 0
                pathway_source_databases = set()
                unit_canonicalization_statuses = []
                mapping_signatures = {}
                for collection_name, genes, matched in (
                    ('matched_genes', matched_genes, True),
                    ('unmatched_genes', unmatched_genes, False),
                ):
                    for index, gene in enumerate(genes):
                        if not isinstance(gene, dict):
                            raise ValueError(
                                f"{collection_name}[{index}] is invalid for "
                                f"{pathway_id}"
                            )
                        validation = self._validate_canonical_gene_unit(
                            gene, pathway_id, matched=matched
                        )
                        canonical_keys.append(validation['gene_key'])
                        source_xref_keys.extend(validation['source_xref_keys'])
                        pathway_source_databases.update(
                            validation['source_databases']
                        )
                        unit_canonicalization_statuses.append(
                            validation['canonicalization_complete']
                        )
                        evaluable_node_count += validation[
                            'node_occurrence_count'
                        ]
                        guide_positive_count += int(
                            matched and validation['has_guides']
                        )
                        if matched:
                            mapping_signatures[validation['gene_key']] = (
                                self._mapping_signature(gene, pathway_id)
                            )

                if (
                    len(canonical_keys) != counts['total_genes']
                    or len(set(canonical_keys)) != len(canonical_keys)
                ):
                    raise ValueError(
                        f"Canonical gene keys do not reconcile for {pathway_id}"
                    )
                if (
                    len(source_xref_keys) != counts['unique_xref_count']
                    or len(set(source_xref_keys)) != len(source_xref_keys)
                ):
                    raise ValueError(
                        f"Source xrefs do not reconcile for {pathway_id}"
                    )
                if evaluable_node_count != counts['evaluable_gene_node_count']:
                    raise ValueError(
                        f"Evaluable node occurrences do not reconcile for "
                        f"{pathway_id}"
                    )
                if guide_positive_count != counts['genes_with_guides']:
                    raise ValueError(
                        f"Guide-positive genes do not reconcile for {pathway_id}"
                    )
                if (
                    len(pathway_source_databases) > 1
                    and not all(unit_canonicalization_statuses)
                ):
                    raise ValueError(
                        f"Cross-database canonicalization does not reconcile "
                        f"for {pathway_id}"
                    )
                expected_mapping_by_pathway[pathway_id] = mapping_signatures
                expected_pathway_names[pathway_id] = row.get(
                    'pathway_name', ''
                )

                try:
                    pw_kfi_score = float(row.get('pw_kfi_score', ''))
                    coverage_rate = float(row.get('coverage_rate', ''))
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        f"Invalid pathway rate for {pathway_id}"
                    ) from exc
                if not math.isfinite(pw_kfi_score) or not math.isfinite(
                    coverage_rate
                ):
                    raise ValueError(f"Non-finite pathway rate for {pathway_id}")
                expected_kfi = (
                    counts['genes_with_guides'] / counts['total_genes']
                    if counts['total_genes'] else 0.0
                )
                expected_coverage = (
                    counts['matched_gene_count'] / counts['total_genes']
                    if counts['total_genes'] else 0.0
                )
                if abs(pw_kfi_score - expected_kfi) > 1e-12:
                    raise ValueError(
                        f"PW-KFI does not reconcile for {pathway_id}"
                    )
                if abs(coverage_rate - expected_coverage) > 1e-12:
                    raise ValueError(
                        f"Annotation match coverage does not reconcile for {pathway_id}"
                    )
                pathways.append(row)

        mapping_path = os.path.join(self.intermediate_dir, 'pathway_gene_mapping.tsv')
        pathway_genes = {}
        if not os.path.exists(mapping_path):
            raise FileNotFoundError(
                f"Canonical pathway gene mapping not found.: {mapping_path}"
            )
        required_mapping_columns = {
            'pathway_id', 'pathway_name', 'pathway_gene_database',
            'pathway_gene_id', 'pathway_gene_name', 'gene_key',
            'pw_kfi_gene_key', 'pathway_xrefs', 'source_databases',
            'node_occurrence_count', 'unique_xref_count',
            'canonical_subject_id', 'canonicalization_methods',
            'canonicalization_complete', 'annotation_source_database',
            'matched_annotation_key', 'matched_annotation_id',
            'resolved_annotation_subject_id', 'match_method',
            'mapping_source', 'sgd_id', 'ensembl_id', 'gene_sgd_id',
            'unique_id', 'gene_name', 'has_guides', 'guide_count',
        }
        with open(mapping_path, 'r', encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            missing_columns = required_mapping_columns - set(
                reader.fieldnames or []
            )
            if missing_columns:
                raise ValueError(
                    "Legacy or incomplete pathway mapping cannot be reported. "
                    "Regenerate pathway analysis; missing columns: "
                    + ', '.join(sorted(missing_columns))
                )
            for row in reader:
                pathway_id = row['pathway_id']
                if pathway_id not in expected_mapping_by_pathway:
                    raise ValueError(
                        f"Mapping references an unknown pathway: {pathway_id}"
                    )
                if row.get('pathway_name', '') != expected_pathway_names[pathway_id]:
                    raise ValueError(
                        f"Pathway mapping name does not reconcile for {pathway_id}"
                    )
                gene_key = row.get('pw_kfi_gene_key', '')
                expected_signature = expected_mapping_by_pathway[pathway_id].get(
                    gene_key
                )
                if expected_signature is None:
                    raise ValueError(
                        f"Pathway mapping key does not reconcile for "
                        f"{pathway_id}: {gene_key!r}"
                    )
                if self._mapping_signature(row, pathway_id) != expected_signature:
                    raise ValueError(
                        f"Pathway mapping payload does not reconcile for "
                        f"{pathway_id}: {gene_key}"
                    )
                pathway_genes.setdefault(pathway_id, []).append(row)

        for pathway_id, expected_mapping in expected_mapping_by_pathway.items():
            mapping_keys = [
                row.get('pw_kfi_gene_key', '')
                for row in pathway_genes.get(pathway_id, [])
            ]
            if (
                len(mapping_keys) != len(set(mapping_keys))
                or set(mapping_keys) != set(expected_mapping)
            ):
                raise ValueError(
                    f"Pathway mapping does not reconcile for {pathway_id}"
                )
        
        self.logger.info(f"Pathway data loading complete.: {len(pathways)} Pathway, {len(pathway_genes)} associated genes")
        return pathways, pathway_genes

    def load_pathway_gene_match_audit(self):
        """Load exact-match and non-evaluable pathway gene audit records."""
        audit_path = os.path.join(
            self.intermediate_dir, 'pathway_gene_match_audit.tsv'
        )
        if not os.path.exists(audit_path):
            self.logger.warning(
                f"Pathway gene match audit not found; legacy report data only: "
                f"{audit_path}"
            )
            return []

        required_columns = {
            'pathway_id', 'pathway_gene_database', 'pathway_gene_id',
            'pw_kfi_gene_key', 'node_occurrence_count',
            'denominator_contribution', 'included_in_denominator',
            'status', 'reason', 'canonicalization_complete',
            'pathway_canonicalization_status',
        }
        with open(audit_path, 'r', encoding='utf-8', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            missing_columns = required_columns - set(reader.fieldnames or [])
            if missing_columns:
                raise ValueError(
                    "Pathway gene match audit is missing required columns: "
                    + ', '.join(sorted(missing_columns))
                )
            return list(reader)
    
    def load_annotation_and_guide_data(self):
        """Load legacy annotation data with deterministic conflict handling."""
        annotation_file = os.path.join(self.intermediate_dir, 'annotation_classified.txt')
        
        if not os.path.exists(annotation_file):
            raise FileNotFoundError(f"Annotation result file not found.: {annotation_file}")
        
        records_by_sgd = {}
        sgd_ids_by_unique_id = {}
        
        with open(annotation_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                sgd_id = str(row.get('subject_id', '') or '').strip()
                unique_id = str(row.get('unique_id', '') or '').strip()
                clean_sgd_id = (
                    sgd_id[len('SGDID:'):] if sgd_id.startswith('SGDID:')
                    else sgd_id
                )
                if not clean_sgd_id:
                    continue
                records_by_sgd.setdefault(clean_sgd_id, []).append(row)
                if unique_id:
                    sgd_ids_by_unique_id.setdefault(unique_id, set()).add(
                        clean_sgd_id
                    )

        annotation_data = {}
        guide_data = {}
        for clean_sgd_id in sorted(records_by_sgd):
            records = records_by_sgd[clean_sgd_id]
            identities = sorted({
                str(record.get('unique_id', '') or '').strip()
                for record in records
                if str(record.get('unique_id', '') or '').strip()
            })
            if len(identities) > 1:
                raise ValueError(
                    "Ambiguous legacy annotation subject mapping: "
                    f"{clean_sgd_id} -> {identities}"
                )
            representative = min(
                records,
                key=lambda record: tuple(
                    (key, str(record.get(key, '') or ''))
                    for key in sorted(record)
                ),
            ).copy()
            annotation_data[clean_sgd_id] = representative

            guide_entries = sorted({
                entry.strip()
                for record in records
                for entry in str(
                    record.get('guide_targets', '') or ''
                ).split(',')
                if entry.strip()
            })
            guides = []
            for guide_entry in guide_entries:
                parts = guide_entry.split('/')
                if len(parts) >= 6:
                    guides.append({
                        'unique_id': parts[0],
                        'start': parts[1],
                        'end': parts[2],
                        'sequence': parts[3],
                        'MM0': parts[4],
                        'MM1-2': parts[5],
                        'off_target_match_score': (
                            parts[6] if len(parts) > 6 else ''
                        ),
                        'strand': parts[7] if len(parts) > 7 else '',
                    })
            guide_data[clean_sgd_id] = guides

        unique_id_to_sgd = {}
        for unique_id in sorted(sgd_ids_by_unique_id):
            sgd_ids = sorted(sgd_ids_by_unique_id[unique_id])
            if len(sgd_ids) != 1:
                raise ValueError(
                    "Ambiguous legacy unique_id mapping: "
                    f"{unique_id} -> {sgd_ids}"
                )
            unique_id_to_sgd[unique_id] = sgd_ids[0]
        
        self.logger.info(f"Annotation guide data loading complete.: {len(annotation_data)} Gene, {len(unique_id_to_sgd)} unique_id Mapping")
        return annotation_data, guide_data, unique_id_to_sgd
    
    def create_pathway_database_tables(self):
        """Create a pathway table in the SQLite database."""
        db_path = os.path.join(self.web_data_dir, 'twheel_data.db')
        
        with sqlite3.connect(db_path) as conn:
            cursor = conn.cursor()
            
            # パスウェイ情報テーブル作成
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS pathway_info (
                    id INTEGER PRIMARY KEY,
                    pathway_id TEXT UNIQUE,
                    pathway_name TEXT,
                    pathway_description TEXT,
                    pathway_organism TEXT,
                    total_gene_nodes INTEGER,
                    evaluable_gene_nodes INTEGER,
                    unique_xref_count INTEGER,
                    total_genes INTEGER,
                    gene_count_basis TEXT,
                    matched_genes INTEGER,
                    genes_with_guides INTEGER,
                    pw_kfi_score REAL,
                    targetable_type TEXT,
                    coverage_rate REAL,
                    canonicalization_status TEXT,
                    gpml_file_path TEXT
                )
            ''')
            
            # パスウェイ-遺伝子関係テーブル作成
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS pathway_gene_relations (
                    id INTEGER PRIMARY KEY,
                    pathway_id TEXT,
                    pathway_gene_database TEXT,
                    pathway_gene_id TEXT,
                    gene_key TEXT,
                    pw_kfi_gene_key TEXT,
                    pathway_xrefs TEXT,
                    source_databases TEXT,
                    node_occurrence_count INTEGER,
                    unique_xref_count INTEGER,
                    canonical_subject_id TEXT,
                    canonicalization_methods TEXT,
                    canonicalization_complete BOOLEAN,
                    annotation_source_database TEXT,
                    matched_annotation_key TEXT,
                    matched_annotation_id TEXT,
                    resolved_annotation_subject_id TEXT,
                    match_method TEXT,
                    mapping_source TEXT,
                    sgd_id TEXT,
                    ensembl_id TEXT,
                    unique_id TEXT,
                    gene_sgd_id TEXT,
                    gene_name TEXT,
                    has_guides BOOLEAN,
                    guide_count INTEGER,
                    FOREIGN KEY (pathway_id) REFERENCES pathway_info (pathway_id)
                )
            ''')

            cursor.execute('''
                CREATE TABLE IF NOT EXISTS pathway_gene_match_audit (
                    id INTEGER PRIMARY KEY,
                    pathway_id TEXT,
                    pathway_name TEXT,
                    gpml_file_path TEXT,
                    pathway_gene_database TEXT,
                    pathway_gene_id TEXT,
                    pathway_gene_name TEXT,
                    gene_key TEXT,
                    pw_kfi_gene_key TEXT,
                    annotation_key TEXT,
                    node_occurrence_count INTEGER,
                    denominator_contribution INTEGER,
                    included_in_denominator BOOLEAN,
                    status TEXT,
                    reason TEXT,
                    canonical_subject_id TEXT,
                    canonicalization_methods TEXT,
                    canonicalization_complete BOOLEAN,
                    pathway_canonicalization_status TEXT,
                    match_method TEXT,
                    mapping_source TEXT,
                    annotation_source_database TEXT,
                    matched_annotation_key TEXT,
                    matched_annotation_id TEXT,
                    resolved_annotation_subject_id TEXT,
                    sgd_id TEXT,
                    ensembl_id TEXT,
                    unique_id TEXT
                )
            ''')
            cursor.execute('''
                CREATE INDEX IF NOT EXISTS idx_pathway_gene_match_audit_pathway
                ON pathway_gene_match_audit (pathway_id)
            ''')
            cursor.execute('''
                CREATE INDEX IF NOT EXISTS idx_pathway_gene_match_audit_status
                ON pathway_gene_match_audit (status, reason)
            ''')

            # Migrate an existing report database without discarding user data.
            required_relation_columns = {
                'pathway_gene_database': 'TEXT',
                'pathway_gene_id': 'TEXT',
                'gene_key': 'TEXT',
                'pw_kfi_gene_key': 'TEXT',
                'pathway_xrefs': 'TEXT',
                'source_databases': 'TEXT',
                'node_occurrence_count': 'INTEGER',
                'unique_xref_count': 'INTEGER',
                'canonical_subject_id': 'TEXT',
                'canonicalization_methods': 'TEXT',
                'canonicalization_complete': 'BOOLEAN',
                'annotation_source_database': 'TEXT',
                'matched_annotation_key': 'TEXT',
                'matched_annotation_id': 'TEXT',
                'resolved_annotation_subject_id': 'TEXT',
                'match_method': 'TEXT',
                'mapping_source': 'TEXT',
                'sgd_id': 'TEXT',
                'ensembl_id': 'TEXT',
                'unique_id': 'TEXT',
            }
            required_info_columns = {
                'total_gene_nodes': 'INTEGER',
                'evaluable_gene_nodes': 'INTEGER',
                'unique_xref_count': 'INTEGER',
                'gene_count_basis': 'TEXT',
                'canonicalization_status': 'TEXT',
            }
            required_audit_columns = {
                'pw_kfi_gene_key': 'TEXT',
                'node_occurrence_count': 'INTEGER',
                'denominator_contribution': 'INTEGER',
                'canonical_subject_id': 'TEXT',
                'canonicalization_methods': 'TEXT',
                'canonicalization_complete': 'BOOLEAN',
                'pathway_canonicalization_status': 'TEXT',
            }
            for table_name, required_columns in (
                ('pathway_info', required_info_columns),
                ('pathway_gene_relations', required_relation_columns),
                ('pathway_gene_match_audit', required_audit_columns),
            ):
                existing_columns = {
                    row[1]
                    for row in cursor.execute(
                        f'PRAGMA table_info({table_name})'
                    ).fetchall()
                }
                for column_name, column_type in required_columns.items():
                    if column_name not in existing_columns:
                        cursor.execute(
                            f'ALTER TABLE {table_name} '
                            f'ADD COLUMN {column_name} {column_type}'
                        )
            
            conn.commit()
        
        self.logger.info(f"Pathway database table creation complete.: {db_path}")
        return db_path
    
    def populate_pathway_database(
        self, pathways, pathway_genes, db_path, pathway_gene_audit=None
    ):
        """Insert pathway data into the database."""
        with sqlite3.connect(db_path) as conn:
            cursor = conn.cursor()
            
            # 既存のパスウェイデータをクリア
            cursor.execute('DELETE FROM pathway_gene_relations')
            cursor.execute('DELETE FROM pathway_gene_match_audit')
            cursor.execute('DELETE FROM pathway_info')
            
            # パスウェイ情報を挿入
            matched_counts_by_pathway = {}
            for pathway in pathways:
                matched_gene_count_value = pathway.get('matched_gene_count', '')
                if str(matched_gene_count_value).strip():
                    try:
                        matched_gene_count = int(matched_gene_count_value)
                    except (TypeError, ValueError) as exc:
                        raise ValueError(
                            f"Invalid matched_gene_count for pathway "
                            f"{pathway.get('pathway_id', '')}: "
                            f"{matched_gene_count_value!r}"
                        ) from exc
                else:
                    if 'matched_genes' not in pathway:
                        raise ValueError(
                            f"Missing matched gene count and list for pathway "
                            f"{pathway.get('pathway_id', '')}"
                        )
                    try:
                        parsed_matched_genes = self._parse_matched_genes(
                            pathway['matched_genes']
                        )
                    except (
                        TypeError, ValueError, SyntaxError, json.JSONDecodeError
                    ) as exc:
                        raise ValueError(
                            f"Cannot parse matched genes for pathway "
                            f"{pathway.get('pathway_id', '')}"
                        ) from exc
                    if not isinstance(parsed_matched_genes, list):
                        raise ValueError(
                            f"Matched genes must be a list for pathway "
                            f"{pathway.get('pathway_id', '')}"
                        )
                    matched_gene_count = len(parsed_matched_genes)

                if matched_gene_count < 0:
                    raise ValueError(
                        f"Negative matched_gene_count for pathway "
                        f"{pathway.get('pathway_id', '')}"
                    )
                matched_counts_by_pathway[pathway['pathway_id']] = (
                    matched_gene_count
                )
                
                cursor.execute('''
                    INSERT INTO pathway_info (
                        pathway_id, pathway_name, pathway_description, pathway_organism,
                        total_gene_nodes, evaluable_gene_nodes, unique_xref_count,
                        total_genes, gene_count_basis, matched_genes,
                        genes_with_guides, pw_kfi_score,
                        targetable_type, coverage_rate, canonicalization_status,
                        gpml_file_path
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    pathway['pathway_id'],
                    pathway['pathway_name'],
                    '',  # pathway_description は後で追加可能
                    pathway.get('pathway_organism', ''),
                    int(pathway.get('total_gene_nodes', 0) or 0),
                    int(pathway.get('evaluable_gene_node_count', 0) or 0),
                    int(pathway.get('unique_xref_count', 0) or 0),
                    int(pathway.get('total_genes', 0)),
                    pathway.get('gene_count_basis', ''),
                    matched_gene_count,
                    int(pathway.get('genes_with_guides', 0)),
                    float(pathway.get('pw_kfi_score', 0)),
                    pathway.get('targetable_type', ''),
                    float(pathway.get('coverage_rate', 0)),
                    pathway.get('canonicalization_status', ''),
                    pathway.get('gpml_file_path', '')
                ))

            unknown_mapping_pathways = (
                set(pathway_genes) - set(matched_counts_by_pathway)
            )
            if unknown_mapping_pathways:
                raise ValueError(
                    "Pathway mapping contains unknown pathways: "
                    + ', '.join(sorted(unknown_mapping_pathways))
                )
            for pathway_id, matched_gene_count in (
                matched_counts_by_pathway.items()
            ):
                mapping_keys = [
                    str(gene.get('pw_kfi_gene_key', '') or '')
                    for gene in pathway_genes.get(pathway_id, [])
                ]
                if (
                    any(not key for key in mapping_keys)
                    or len(mapping_keys) != len(set(mapping_keys))
                    or len(mapping_keys) != matched_gene_count
                ):
                    raise ValueError(
                        "Canonical pathway mapping does not match "
                        f"matched_gene_count for {pathway_id}"
                    )
            
            # パスウェイ-遺伝子関係を挿入
            for pathway_id, genes in pathway_genes.items():
                for gene in genes:
                    # guide_count の安全な変換
                    guide_count = gene.get('guide_count', 0)
                    if isinstance(guide_count, str):
                        try:
                            guide_count = int(guide_count)
                        except ValueError:
                            guide_count = 0
                    elif not isinstance(guide_count, int):
                        guide_count = 0
                    
                    cursor.execute('''
                        INSERT INTO pathway_gene_relations (
                            pathway_id, pathway_gene_database, pathway_gene_id,
                            gene_key, pw_kfi_gene_key, pathway_xrefs,
                            source_databases, node_occurrence_count,
                            unique_xref_count, canonical_subject_id,
                            canonicalization_methods,
                            canonicalization_complete,
                            annotation_source_database,
                            matched_annotation_key, matched_annotation_id,
                            resolved_annotation_subject_id, match_method,
                            mapping_source, sgd_id, ensembl_id, unique_id,
                            gene_sgd_id, gene_name, has_guides, guide_count
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ''', (
                        pathway_id,
                        gene.get('pathway_gene_database', ''),
                        gene.get('pathway_gene_id', ''),
                        gene.get('gene_key', ''),
                        gene.get('pw_kfi_gene_key', ''),
                        gene.get('pathway_xrefs', ''),
                        gene.get('source_databases', ''),
                        int(gene.get('node_occurrence_count', 0) or 0),
                        int(gene.get('unique_xref_count', 0) or 0),
                        gene.get('canonical_subject_id', ''),
                        gene.get('canonicalization_methods', ''),
                        str(gene.get(
                            'canonicalization_complete', 'False'
                        )).lower() == 'true',
                        gene.get('annotation_source_database', ''),
                        gene.get('matched_annotation_key', ''),
                        gene.get('matched_annotation_id', ''),
                        gene.get('resolved_annotation_subject_id', ''),
                        gene.get('match_method', ''),
                        gene.get('mapping_source', ''),
                        gene.get('sgd_id', ''),
                        gene.get('ensembl_id', ''),
                        gene.get('unique_id', ''),
                        gene.get('gene_sgd_id', ''),
                        gene.get('gene_name', ''),
                        str(gene.get('has_guides', 'False')).lower() == 'true',
                        guide_count
                    ))

            for audit_row in pathway_gene_audit or []:
                cursor.execute('''
                    INSERT INTO pathway_gene_match_audit (
                        pathway_id, pathway_name, gpml_file_path,
                        pathway_gene_database, pathway_gene_id,
                        pathway_gene_name, gene_key, annotation_key,
                        pw_kfi_gene_key, node_occurrence_count,
                        denominator_contribution, included_in_denominator,
                        status, reason, canonical_subject_id,
                        canonicalization_methods, canonicalization_complete,
                        pathway_canonicalization_status, match_method,
                        mapping_source, annotation_source_database,
                        matched_annotation_key, matched_annotation_id,
                        resolved_annotation_subject_id, sgd_id, ensembl_id,
                        unique_id
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    audit_row.get('pathway_id', ''),
                    audit_row.get('pathway_name', ''),
                    audit_row.get('gpml_file_path', ''),
                    audit_row.get('pathway_gene_database', ''),
                    audit_row.get('pathway_gene_id', ''),
                    audit_row.get('pathway_gene_name', ''),
                    audit_row.get('gene_key', ''),
                    audit_row.get('annotation_key', ''),
                    audit_row.get('pw_kfi_gene_key', ''),
                    int(audit_row.get('node_occurrence_count', 0) or 0),
                    int(audit_row.get('denominator_contribution', 0) or 0),
                    str(audit_row.get(
                        'included_in_denominator', 'False'
                    )).lower() == 'true',
                    audit_row.get('status', ''),
                    audit_row.get('reason', ''),
                    audit_row.get('canonical_subject_id', ''),
                    audit_row.get('canonicalization_methods', ''),
                    str(audit_row.get(
                        'canonicalization_complete', 'False'
                    )).lower() == 'true',
                    audit_row.get('pathway_canonicalization_status', ''),
                    audit_row.get('match_method', ''),
                    audit_row.get('mapping_source', ''),
                    audit_row.get('annotation_source_database', ''),
                    audit_row.get('matched_annotation_key', ''),
                    audit_row.get('matched_annotation_id', ''),
                    audit_row.get('resolved_annotation_subject_id', ''),
                    audit_row.get('sgd_id', ''),
                    audit_row.get('ensembl_id', ''),
                    audit_row.get('unique_id', ''),
                ))
            
            conn.commit()
        
        self.logger.info("Pathway data insertion into the database is complete.")
    
    def parse_gpml_file(self, gpml_file_path):
        """
        Analyze GPML files to extract pathway diagram data.
        
        Args:
            gpml_file_path (str): Path to the GPML file
            
        Returns:
            dict: Pathway diagram data (nodes, edges, metadata)
        """
        try:
            if not os.path.exists(gpml_file_path):
                self.logger.warning(f"GPML file not found.: {gpml_file_path}")
                return None
            
            # XMLファイルを解析
            tree = ET.parse(gpml_file_path)
            root = tree.getroot()
            
            # 名前空間を取得
            ns = {'gpml': 'http://pathvisio.org/GPML/2013a'} if root.tag.startswith('{') else {}
            
            # パスウェイメタデータを抽出
            pathway_name = root.get('Name', 'Unknown Pathway')
            organism = root.get('Organism', 'Unknown Organism')
            
            # Graphics要素からキャンバスサイズを取得
            graphics_elem = root.find('.//gpml:Graphics', ns) if ns else root.find('.//Graphics')
            board_width = float(graphics_elem.get('BoardWidth', 800)) if graphics_elem is not None else 800
            board_height = float(graphics_elem.get('BoardHeight', 600)) if graphics_elem is not None else 600
            
            # DataNodeを抽出（遺伝子・代謝物ノード）
            nodes = []
            data_nodes = root.findall('.//gpml:DataNode', ns) if ns else root.findall('.//DataNode')
            
            for node in data_nodes:
                node_type = node.get('Type', 'Unknown')
                node_id = node.get('GraphId', '')
                text_label = node.get('TextLabel', '')
                
                # Graphics情報を取得
                graphics = node.find('gpml:Graphics', ns) if ns else node.find('Graphics')
                center_x = float(graphics.get('CenterX', 0)) if graphics is not None else 0
                center_y = float(graphics.get('CenterY', 0)) if graphics is not None else 0
                width = float(graphics.get('Width', 80)) if graphics is not None else 80
                height = float(graphics.get('Height', 20)) if graphics is not None else 20
                
                # Xref情報を取得（SGD IDなど）
                xref = node.find('gpml:Xref', ns) if ns else node.find('Xref')
                database = (
                    xref.get('Database', '').strip() if xref is not None else ''
                )
                xref_id = xref.get('ID', '').strip() if xref is not None else ''
                
                node_data = {
                    'id': node_id,
                    'type': node_type,
                    'label': text_label,
                    'x': center_x,
                    'y': center_y,
                    'width': width,
                    'height': height,
                    'database': database,
                    'xref_id': xref_id
                }
                
                nodes.append(node_data)
            
            # Interactionを抽出（相互作用エッジ）
            edges = []
            interactions = root.findall('.//gpml:Interaction', ns) if ns else root.findall('.//Interaction')
            
            for interaction in interactions:
                interaction_id = interaction.get('GraphId', '')
                
                # Graphics情報を取得
                graphics = interaction.find('gpml:Graphics', ns) if ns else interaction.find('Graphics')
                connector_type = graphics.get('ConnectorType', 'Straight') if graphics is not None else 'Straight'
                
                # Point情報を取得
                points = []
                if graphics is not None:
                    point_elements = graphics.findall('gpml:Point', ns) if ns else graphics.findall('Point')
                    for point in point_elements:
                        x = float(point.get('X', 0))
                        y = float(point.get('Y', 0))
                        graph_ref = point.get('GraphRef', '')
                        points.append({'x': x, 'y': y, 'ref': graph_ref})
                
                edge_data = {
                    'id': interaction_id,
                    'type': 'interaction',
                    'connector_type': connector_type,
                    'points': points
                }
                
                edges.append(edge_data)
            
            pathway_data = {
                'pathway_name': pathway_name,
                'organism': organism,
                'board_width': board_width,
                'board_height': board_height,
                'nodes': nodes,
                'edges': edges,
                'gpml_file': gpml_file_path
            }
            
            self.logger.debug(f"GPML analysis complete: {pathway_name} - {len(nodes)} node, {len(edges)} edge")
            return pathway_data
            
        except Exception as e:
            self.logger.error(f"GPML parsing error ({gpml_file_path}): {e}")
            return None
    
    def generate_pathway_svg(self, pathway_data, guide_status_map=None):
        """
        Generate SVG strings from pathway data.
        
        Args:
            pathway_data (dict): Data returned from parse_gpml_file
            guide_status_map (dict): Exact ``Database::ID`` -> guide status mapping.
            
        Returns:
            str: SVG string
        """
        if not pathway_data:
            return None
        
        board_width = pathway_data['board_width']
        board_height = pathway_data['board_height']
        nodes = pathway_data['nodes']
        edges = pathway_data['edges']
        
        # SVGヘッダー
        svg_parts = [
            f'<svg width="{board_width}" height="{board_height}" viewBox="0 0 {board_width} {board_height}" xmlns="http://www.w3.org/2000/svg">',
            '<defs>',
            '<style>',
            '.gene-node { stroke-width: 2; cursor: pointer; }',
            '.gene-with-guides { fill: #28a745; stroke: #1e7e34; }',
            '.gene-without-guides { fill: #dc3545; stroke: #bd2130; }',
            '.gene-unknown { fill: #6c757d; stroke: #545b62; }',
            '.metabolite-node { fill: #17a2b8; stroke: #117a8b; stroke-width: 2; }',
            '.interaction-line { stroke: #333; stroke-width: 1.5; fill: none; }',
            '.node-label { font-family: Arial, sans-serif; font-size: 10px; text-anchor: middle; pointer-events: none; }',
            '</style>',
            '</defs>'
        ]
        
        # エッジ（相互作用線）を先に描画
        for edge in edges:
            points = edge['points']
            if len(points) >= 2:
                path_data = f"M {points[0]['x']} {points[0]['y']}"
                for point in points[1:]:
                    path_data += f" L {point['x']} {point['y']}"
                
                svg_parts.append(f'<path d="{path_data}" class="interaction-line" />')
        
        # ノード（遺伝子・代謝物）を描画
        for node in nodes:
            x, y = node['x'], node['y']
            width, height = node['width'], node['height']
            node_type = node['type']
            label = node['label']
            xref_id = node['xref_id']
            database = node['database']
            gene_key = f"{database}::{xref_id}"
            escaped_label = html.escape(str(label), quote=True)
            escaped_database = html.escape(str(database), quote=True)
            escaped_xref_id = html.escape(str(xref_id), quote=True)
            escaped_gene_key = html.escape(gene_key, quote=True)
            
            # Color only an exact database-and-identifier match.
            if node_type == 'GeneProduct':
                gene_status = (
                    guide_status_map.get(gene_key)
                    if guide_status_map else None
                )
                if (
                    gene_status is not None
                    and gene_status.get('match_status', 'matched') == 'matched'
                ):
                    if gene_status['has_guides']:
                        node_class = 'gene-with-guides'
                    else:
                        node_class = 'gene-without-guides'
                else:
                    node_class = 'gene-unknown'
                
                # 遺伝子ノード（楕円形）
                svg_parts.append(
                    f'<ellipse cx="{x}" cy="{y}" rx="{width/2}" ry="{height/2}" '
                    f'class="gene-node {node_class}" '
                    f'data-gene-database="{escaped_database}" '
                    f'data-gene-id="{escaped_xref_id}" '
                    f'data-gene-key="{escaped_gene_key}" '
                    f'data-sgd-id="{html.escape(str((gene_status or {}).get("sgd_id", "")), quote=True)}" '
                    f'data-gene-name="{escaped_label}" '
                    f'title="{escaped_label} ({escaped_database}:{escaped_xref_id})" />'
                )
            
            elif node_type == 'Metabolite':
                # 代謝物ノード（四角形）
                svg_parts.append(
                    f'<rect x="{x - width/2}" y="{y - height/2}" width="{width}" height="{height}" '
                    f'class="metabolite-node" title="{label}" />'
                )
            
            else:
                # その他のノード（四角形、グレー）
                svg_parts.append(
                    f'<rect x="{x - width/2}" y="{y - height/2}" width="{width}" height="{height}" '
                    f'fill="#e9ecef" stroke="#adb5bd" stroke-width="1" title="{label}" />'
                )
            
            # ラベルテキスト
            # 長いラベルは切り詰める
            display_label = label[:12] + '...' if len(label) > 15 else label
            svg_parts.append(
                f'<text x="{x}" y="{y + 3}" class="node-label" fill="#000">'
                f'{html.escape(str(display_label))}</text>'
            )
        
        # 凡例を追加
        legend_y = board_height - 80
        svg_parts.extend([
            f'<g id="legend">',
            f'<rect x="10" y="{legend_y}" width="200" height="70" fill="white" stroke="#ccc" stroke-width="1" opacity="0.9" />',
            f'<text x="20" y="{legend_y + 15}" class="node-label" font-weight="bold">Legend</text>',
            f'<ellipse cx="30" cy="{legend_y + 30}" rx="8" ry="6" class="gene-node gene-with-guides" />',
            f'<text x="45" y="{legend_y + 33}" class="node-label" font-size="9px">With Guides</text>',
            f'<ellipse cx="30" cy="{legend_y + 45}" rx="8" ry="6" class="gene-node gene-without-guides" />',
            f'<text x="45" y="{legend_y + 48}" class="node-label" font-size="9px">Without Guides</text>',
            f'<ellipse cx="125" cy="{legend_y + 30}" rx="8" ry="6" class="gene-node gene-unknown" />',
            f'<text x="140" y="{legend_y + 33}" class="node-label" font-size="9px">Unknown</text>',
            f'<rect x="120" y="{legend_y + 39}" width="16" height="12" class="metabolite-node" />',
            f'<text x="140" y="{legend_y + 48}" class="node-label" font-size="9px">Metabolite</text>',
            f'</g>'
        ])
        
        svg_parts.append('</svg>')
        
        return '\n'.join(svg_parts)
    
    @staticmethod
    def _parse_matched_genes(matched_genes_value):
        """Parse the serialized matched-gene list without interpreting free text."""
        if isinstance(matched_genes_value, list):
            return matched_genes_value
        if not isinstance(matched_genes_value, str):
            raise TypeError('Serialized pathway gene data must be a string or list')
        if not matched_genes_value.strip():
            return []
        if matched_genes_value.startswith('[') and matched_genes_value.endswith(']'):
            return ast.literal_eval(matched_genes_value)
        return json.loads(matched_genes_value)

    @staticmethod
    def _exact_pathway_identifier(gene_info):
        """Return an explicit pathway database and ID, with SGD-only legacy support."""
        database = (gene_info.get('pathway_gene_database') or '').strip()
        gene_id = (gene_info.get('pathway_gene_id') or '').strip()
        if database and gene_id:
            return database, gene_id

        # Older pathway output was produced only from an exact SGD GPML ID.
        legacy_id = (gene_info.get('original_gene_id') or '').strip()
        if re.fullmatch(r'S\d{9}', legacy_id):
            return 'SGD', legacy_id
        return '', ''

    @classmethod
    def _pathway_xrefs(cls, gene_info):
        """Return every exact source xref represented by one PW-KFI gene."""
        value = gene_info.get('pathway_xrefs', [])
        if isinstance(value, str) and value.strip():
            value = cls._parse_matched_genes(value)
        if value:
            if not isinstance(value, list):
                raise ValueError('pathway_xrefs must be a list')
            xrefs = []
            for index, xref in enumerate(value):
                if not isinstance(xref, dict):
                    raise ValueError(f'pathway_xrefs[{index}] must be an object')
                database = str(xref.get('database', '') or '').strip()
                gene_id = str(xref.get('gene_id', '') or '').strip()
                if not database or not gene_id:
                    raise ValueError(
                        f'pathway_xrefs[{index}] lacks database or gene_id'
                    )
                expected_gene_key = f'{database}::{gene_id}'
                supplied_gene_key = str(
                    xref.get('gene_key', '') or ''
                ).strip()
                if supplied_gene_key and supplied_gene_key != expected_gene_key:
                    raise ValueError(
                        f'pathway_xrefs[{index}] has an inconsistent gene_key'
                    )
                normalized = xref.copy()
                normalized.update({
                    'database': database,
                    'gene_id': gene_id,
                    'gene_key': expected_gene_key,
                })
                xrefs.append(normalized)
            return sorted(
                xrefs, key=lambda xref: (xref['database'], xref['gene_id'])
            )

        database = str(
            gene_info.get('pathway_gene_database')
            or gene_info.get('database')
            or ''
        ).strip()
        gene_id = str(
            gene_info.get('pathway_gene_id')
            or gene_info.get('gene_id')
            or ''
        ).strip()
        if not database or not gene_id:
            database, gene_id = cls._exact_pathway_identifier(gene_info)
        if not database or not gene_id:
            return []
        return [{
            'database': database,
            'gene_id': gene_id,
            'gene_key': f'{database}::{gene_id}',
            'pathway_gene_name': gene_info.get('pathway_gene_name', ''),
            'match_status': gene_info.get('match_status', 'matched'),
            'reason': gene_info.get('reason', ''),
        }]

    @staticmethod
    def _nonnegative_int(value, field_name, pathway_id, default=None):
        """Parse a report count without silently coercing invalid values."""
        if value is None or not str(value).strip():
            if default is not None:
                return default
            raise ValueError(f'Missing {field_name} for {pathway_id}')
        try:
            parsed = int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f'Invalid {field_name} for {pathway_id}: {value!r}'
            ) from exc
        if parsed < 0:
            raise ValueError(f'Negative {field_name} for {pathway_id}')
        return parsed

    @staticmethod
    def _strict_bool(value, field_name, pathway_id):
        """Parse an artifact boolean without treating arbitrary text as false."""
        if isinstance(value, bool):
            return value
        normalized = str(value).strip().lower()
        if normalized == 'true':
            return True
        if normalized == 'false':
            return False
        raise ValueError(
            f'Invalid {field_name} for {pathway_id}: {value!r}'
        )

    @staticmethod
    def _normalized_string_list(value):
        """Normalize list-valued artifact fields written as CSV text."""
        if value is None:
            return ()
        if isinstance(value, str):
            items = value.split(',')
        elif isinstance(value, (list, tuple, set)):
            items = value
        else:
            raise ValueError(
                f'Expected a string or list, got {type(value).__name__}'
            )
        return tuple(sorted({
            str(item).strip() for item in items if str(item).strip()
        }))

    @classmethod
    def _freeze_artifact_value(cls, value):
        """Convert nested provenance data to an order-independent value."""
        if isinstance(value, dict):
            return tuple(sorted(
                (str(key), cls._freeze_artifact_value(item))
                for key, item in value.items()
            ))
        if isinstance(value, (list, tuple, set)):
            frozen = [cls._freeze_artifact_value(item) for item in value]
            return tuple(sorted(frozen, key=repr))
        return value

    @staticmethod
    def _guide_entries(value, pathway_id):
        """Normalize serialized guide entries for artifact reconciliation."""
        if value is None:
            return ()
        if isinstance(value, str):
            items = value.split(',')
        elif isinstance(value, (list, tuple, set)):
            items = value
        else:
            raise ValueError(
                f'Invalid guide_targets for {pathway_id}: {value!r}'
            )
        return tuple(sorted({
            str(item).strip() for item in items if str(item).strip()
        }))

    @classmethod
    def _validate_canonical_gene_unit(cls, gene, pathway_id, matched):
        """Validate one canonical PW-KFI denominator unit and its source nodes."""
        gene_key = str(gene.get('pw_kfi_gene_key', '') or '').strip()
        if not gene_key:
            raise ValueError(
                f'Canonical gene provenance is missing for {pathway_id}'
            )
        xrefs = cls._pathway_xrefs(gene)
        if not xrefs:
            raise ValueError(
                f'Canonical gene source xrefs are missing for {pathway_id}'
            )
        source_xref_keys = [xref['gene_key'] for xref in xrefs]
        if len(source_xref_keys) != len(set(source_xref_keys)):
            raise ValueError(
                f'Duplicate source xrefs in canonical gene {gene_key}'
            )

        node_occurrence_count = sum(
            cls._nonnegative_int(
                xref.get('node_occurrence_count', ''),
                'pathway_xrefs.node_occurrence_count', pathway_id,
            )
            for xref in xrefs
        )
        if node_occurrence_count != cls._nonnegative_int(
            gene.get('node_occurrence_count', ''),
            'node_occurrence_count', pathway_id,
        ):
            raise ValueError(
                f'Node occurrences do not reconcile for {pathway_id}: {gene_key}'
            )
        if len(xrefs) != cls._nonnegative_int(
            gene.get('unique_xref_count', ''),
            'unique_xref_count', pathway_id,
        ):
            raise ValueError(
                f'Unique xrefs do not reconcile for {pathway_id}: {gene_key}'
            )

        primary_xref = xrefs[0]
        expected_primary = (
            primary_xref['database'], primary_xref['gene_id'],
            primary_xref['gene_key'],
        )
        supplied_primary = (
            str(gene.get('pathway_gene_database', '') or '').strip(),
            str(gene.get('pathway_gene_id', '') or '').strip(),
            str(gene.get('gene_key', '') or '').strip(),
        )
        if supplied_primary != expected_primary:
            raise ValueError(
                f'Primary pathway xref does not reconcile for {pathway_id}: '
                f'{gene_key}'
            )
        expected_databases = tuple(sorted({
            xref['database'] for xref in xrefs
        }))
        if cls._normalized_string_list(
            gene.get('source_databases', ())
        ) != expected_databases:
            raise ValueError(
                f'Source databases do not reconcile for {pathway_id}: {gene_key}'
            )
        canonicalization_complete = cls._strict_bool(
            gene.get('canonicalization_complete', ''),
            'canonicalization_complete', pathway_id,
        )

        guides = cls._guide_entries(gene.get('guide_targets', ''), pathway_id)
        has_guides = cls._strict_bool(
            gene.get('has_guides', ''), 'has_guides', pathway_id
        )
        guide_count = cls._nonnegative_int(
            gene.get('guide_count', ''), 'guide_count', pathway_id
        )
        if has_guides != bool(guides) or guide_count != len(guides):
            raise ValueError(
                f'Guide counts do not reconcile for {pathway_id}: {gene_key}'
            )
        if not matched and (has_guides or guide_count or guides):
            raise ValueError(
                f'Unmatched canonical gene has guides for {pathway_id}: {gene_key}'
            )
        return {
            'gene_key': gene_key,
            'source_xref_keys': source_xref_keys,
            'source_databases': expected_databases,
            'node_occurrence_count': node_occurrence_count,
            'canonicalization_complete': canonicalization_complete,
            'has_guides': has_guides,
        }

    @classmethod
    def _mapping_signature(cls, gene, pathway_id):
        """Return the exact mapping payload represented by a matched gene."""
        xrefs = cls._pathway_xrefs(gene)
        return (
            str(gene.get('pathway_gene_database', '') or ''),
            str(gene.get('pathway_gene_id', '') or ''),
            str(gene.get('pathway_gene_name', '') or ''),
            str(gene.get('gene_key', '') or ''),
            str(gene.get('pw_kfi_gene_key', '') or ''),
            cls._freeze_artifact_value(xrefs),
            cls._normalized_string_list(gene.get('source_databases', ())),
            cls._nonnegative_int(
                gene.get('node_occurrence_count', ''),
                'node_occurrence_count', pathway_id,
            ),
            cls._nonnegative_int(
                gene.get('unique_xref_count', ''),
                'unique_xref_count', pathway_id,
            ),
            str(gene.get('canonical_subject_id', '') or ''),
            cls._normalized_string_list(
                gene.get('canonicalization_methods', ())
            ),
            cls._strict_bool(
                gene.get('canonicalization_complete', ''),
                'canonicalization_complete', pathway_id,
            ),
            str(gene.get('annotation_source_database', '') or ''),
            str(gene.get('matched_annotation_key', '') or ''),
            str(gene.get('matched_annotation_id', '') or ''),
            str(gene.get('resolved_annotation_subject_id', '') or ''),
            str(gene.get('match_method', '') or ''),
            str(gene.get('mapping_source', '') or ''),
            str(gene.get('sgd_id', '') or ''),
            str(gene.get('ensembl_id', '') or ''),
            str(
                gene.get('gene_sgd_id', gene.get('sgd_id', '')) or ''
            ),
            str(gene.get('unique_id', '') or ''),
            str(gene.get('gene_name', '') or ''),
            cls._strict_bool(
                gene.get('has_guides', ''), 'has_guides', pathway_id
            ),
            cls._nonnegative_int(
                gene.get('guide_count', ''), 'guide_count', pathway_id
            ),
        )

    def get_pathway_diagram_data(self, pathway_id, pathways, unique_id_to_sgd=None):
        """
        Obtain diagrammatic data for a specific pathway.
        
        Args:
            pathway_id (str): Pathway ID
            pathways (list): Pathway data list
            unique_id_to_sgd (dict): Deprecated legacy mapping; not used for matching.
            
        Returns:
            dict: Pathway diagram data (SVG string, including metadata)
        """
        # 指定されたパスウェイを検索
        target_pathway = None
        for pathway in pathways:
            if pathway['pathway_id'] == pathway_id:
                target_pathway = pathway
                break
        
        if not target_pathway:
            self.logger.warning(f"Pathway not found.: {pathway_id}")
            return None
        
        # GPMLファイルパスを取得
        gpml_file_path = target_pathway.get('gpml_file_path', '')
        if not gpml_file_path or not os.path.exists(gpml_file_path):
            self.logger.warning(f"GPML file not found.: {gpml_file_path}")
            return None
        
        # GPMLファイルを解析
        pathway_data = self.parse_gpml_file(gpml_file_path)
        if not pathway_data:
            return None
        
        # ガイドRNA有無のマッピングを作成
        guide_status_map = {}
        unrenderable_gene_records = []
        matched_genes_value = target_pathway.get('matched_genes', '[]')
        
        try:
            matched_genes = self._parse_matched_genes(matched_genes_value)
            
            if not isinstance(matched_genes, list):
                raise ValueError('matched_genes must be a list')
            for index, gene_info in enumerate(matched_genes):
                if not isinstance(gene_info, dict):
                    raise ValueError(
                        f"matched_genes[{index}] must be an object"
                    )
                xrefs = self._pathway_xrefs(gene_info)
                if not xrefs:
                    raise ValueError(
                        f"matched_genes[{index}] lacks an exact pathway "
                        "database and identifier"
                    )

                unique_id = gene_info.get('unique_id', '')
                guide_targets = str(gene_info.get('guide_targets', '') or '')
                guides = [
                    item for item in guide_targets.split(',') if item.strip()
                ]
                sgd_id = str(gene_info.get('sgd_id', '') or '')
                if sgd_id.startswith('SGDID:'):
                    sgd_id = sgd_id[len('SGDID:'):]
                ensembl_id = str(gene_info.get('ensembl_id', '') or '')
                pw_kfi_gene_key = gene_info.get(
                    'pw_kfi_gene_key', gene_info.get('gene_key', '')
                )
                source_xref_keys = [xref['gene_key'] for xref in xrefs]
                for xref in xrefs:
                    database = xref['database']
                    gene_id = xref['gene_id']
                    gene_key = xref['gene_key']
                    xref_sgd_id = sgd_id or (
                        gene_id if database == 'SGD' else ''
                    )
                    xref_ensembl_id = ensembl_id or (
                        gene_id if database == 'Ensembl' else ''
                    )
                    previous = guide_status_map.get(gene_key)
                    if previous and previous.get(
                        'pw_kfi_gene_key'
                    ) != pw_kfi_gene_key:
                        raise ValueError(
                            f"Source xref {gene_key} maps to multiple PW-KFI genes"
                        )
                    guide_status_map[gene_key] = {
                            'gene_key': gene_key,
                            'pw_kfi_gene_key': pw_kfi_gene_key,
                            'source_xref_keys': source_xref_keys,
                            'pathway_gene_database': database,
                            'pathway_gene_id': gene_id,
                            'annotation_source_database': gene_info.get(
                                'annotation_source_database', 'Unspecified'
                            ),
                            'matched_annotation_key': gene_info.get(
                                'matched_annotation_key', ''
                            ),
                            'matched_annotation_id': gene_info.get(
                                'matched_annotation_id', ''
                            ),
                            'resolved_annotation_subject_id': gene_info.get(
                                'resolved_annotation_subject_id', ''
                            ),
                            'match_method': xref.get(
                                'match_method', gene_info.get('match_method', '')
                            ),
                            'mapping_source': xref.get(
                                'mapping_source', gene_info.get('mapping_source', '')
                            ),
                            'match_status': 'matched',
                            'failure_reason': '',
                            'sgd_id': xref_sgd_id,
                            'ensembl_id': xref_ensembl_id,
                            'has_guides': bool(guides),
                            'guide_count': len(guides),
                            'gene_name': gene_info.get('gene_name', ''),
                            'description': gene_info.get('description', ''),
                            'unique_id': unique_id,
                    }

            for field_name, match_status in (
                ('unmatched_genes', 'unmatched'),
                ('unmatchable_genes', 'not_evaluable'),
            ):
                unavailable_genes = self._parse_matched_genes(
                    target_pathway.get(field_name, '[]')
                )
                if not isinstance(unavailable_genes, list):
                    raise ValueError(f"{field_name} must be a list")
                for index, gene_info in enumerate(unavailable_genes):
                    if not isinstance(gene_info, dict):
                        raise ValueError(
                            f"{field_name}[{index}] must be an object"
                        )
                    xrefs = self._pathway_xrefs(gene_info)
                    if not xrefs:
                        unrenderable_gene_records.append({
                            'source_field': field_name,
                            'index': index,
                            'pathway_gene_database': '',
                            'pathway_gene_id': '',
                            'pathway_gene_name': gene_info.get('gene_name', ''),
                            'match_status': match_status,
                            'failure_reason': gene_info.get(
                                'reason', 'missing_exact_database_or_identifier'
                            ),
                        })
                        continue
                    for xref in xrefs:
                        database = xref['database']
                        gene_id = xref['gene_id']
                        gene_key = xref['gene_key']
                        if gene_key in guide_status_map:
                            continue
                        guide_status_map[gene_key] = {
                        'gene_key': gene_key,
                        'pw_kfi_gene_key': gene_info.get(
                            'pw_kfi_gene_key', gene_key
                        ),
                        'source_xref_keys': [
                            source_xref['gene_key'] for source_xref in xrefs
                        ],
                        'pathway_gene_database': database,
                        'pathway_gene_id': gene_id,
                        'annotation_source_database': '',
                        'matched_annotation_key': '',
                        'matched_annotation_id': '',
                        'resolved_annotation_subject_id': '',
                        'match_method': '',
                        'mapping_source': '',
                        'match_status': match_status,
                        'failure_reason': xref.get(
                            'reason', gene_info.get('reason', match_status)
                        ),
                        'sgd_id': '',
                        'ensembl_id': gene_id if database == 'Ensembl' else '',
                        'has_guides': False,
                        'guide_count': 0,
                        'gene_name': gene_info.get('gene_name', ''),
                        'description': '',
                        'unique_id': '',
                        }

            node_occurrence_counts = Counter(
                f"{node['database']}::{node['xref_id']}"
                for node in pathway_data['nodes']
                if node['type'] == 'GeneProduct'
                and node['database']
                and node['xref_id']
            )
            for gene_key, gene_status in guide_status_map.items():
                gene_status['node_occurrence_count'] = node_occurrence_counts.get(
                    gene_key, 0
                )
        
        except Exception as exc:
            raise ValueError(
                f"Cannot deserialize exact pathway gene data for {pathway_id}: {exc}"
            ) from exc
        
        # SVGを生成
        svg_content = self.generate_pathway_svg(pathway_data, guide_status_map)
        
        if svg_content:
            total_gene_nodes = len([
                node for node in pathway_data['nodes']
                if node['type'] == 'GeneProduct'
            ])
            evaluable_gene_nodes = self._nonnegative_int(
                target_pathway.get('evaluable_gene_node_count', ''),
                'evaluable_gene_node_count', pathway_id,
                default=sum(
                    value.get('node_occurrence_count', 0)
                    for value in guide_status_map.values()
                    if value.get('match_status') in {'matched', 'unmatched'}
                ),
            )
            total_genes = self._nonnegative_int(
                target_pathway.get('total_genes', ''),
                'total_genes', pathway_id,
                default=len({
                    value.get('pw_kfi_gene_key', key)
                    for key, value in guide_status_map.items()
                    if value.get('match_status') in {'matched', 'unmatched'}
                }),
            )
            unique_xref_count = self._nonnegative_int(
                target_pathway.get('unique_xref_count', ''),
                'unique_xref_count', pathway_id,
                default=sum(
                    value.get('match_status') in {'matched', 'unmatched'}
                    for value in guide_status_map.values()
                ),
            )
            genes_with_guides = self._nonnegative_int(
                target_pathway.get('genes_with_guides', ''),
                'genes_with_guides', pathway_id,
                default=len({
                    value.get('pw_kfi_gene_key', key)
                    for key, value in guide_status_map.items()
                    if value.get('match_status') == 'matched'
                    and value.get('has_guides')
                }),
            )
            gene_nodes_with_guides = sum(
                value.get('node_occurrence_count', 0)
                for value in guide_status_map.values()
                if value.get('match_status') == 'matched'
                and value.get('has_guides')
            )

            return {
                'pathway_id': pathway_id,
                'pathway_name': pathway_data['pathway_name'],
                'organism': pathway_data['organism'],
                'svg_content': svg_content,
                'board_width': pathway_data['board_width'],
                'board_height': pathway_data['board_height'],
                'guide_status_map': guide_status_map,
                'unrenderable_gene_records': unrenderable_gene_records,
                'total_genes': total_genes,
                'evaluable_gene_nodes': evaluable_gene_nodes,
                'unique_xref_count': unique_xref_count,
                'total_gene_nodes': total_gene_nodes,
                'genes_with_guides': genes_with_guides,
                'gene_nodes_with_guides': gene_nodes_with_guides,
            }
        
        return None

    def generate_pathway_guide_csv_data(
        self, pathways, unique_id_to_sgd=None, annotation_data=None, guide_data=None
    ):
        """Generate CSV data for pathway-specific guide designs."""
        pathway_guide_data = {}
        
        for pathway in pathways:
            pathway_id = pathway['pathway_id']
            pathway_name = pathway['pathway_name']
            matched_genes_value = pathway.get('matched_genes', '[]')
            
            # matched_genes列をJSONとして解析
            try:
                matched_genes = self._parse_matched_genes(matched_genes_value)
            except Exception as exc:
                raise ValueError(
                    f"Cannot deserialize matched genes for pathway {pathway_id}: "
                    f"{exc}"
                ) from exc
            if not isinstance(matched_genes, list):
                raise ValueError(
                    f"Matched genes must be a list for pathway {pathway_id}"
                )
            
            guides_data = []
            seen_guides = set()
            
            for gene_info in matched_genes:
                if not isinstance(gene_info, dict):
                    raise ValueError(
                        f"Matched gene entries must be objects for pathway "
                        f"{pathway_id}"
                    )
                unique_id = gene_info.get('unique_id', '')
                gene_name = gene_info.get('gene_name', '')
                description = gene_info.get('description', '')
                xrefs = self._pathway_xrefs(gene_info)
                if not xrefs:
                    raise ValueError(
                        f"Pathway {pathway_id} matched gene lacks an exact "
                        "pathway database and identifier"
                    )
                database = xrefs[0]['database']
                gene_id = xrefs[0]['gene_id']
                gene_key = xrefs[0]['gene_key']
                pw_kfi_gene_key = gene_info.get(
                    'pw_kfi_gene_key', gene_key
                )
                source_xref_keys = [xref['gene_key'] for xref in xrefs]
                sgd_id = str(gene_info.get('sgd_id', '') or '')
                if sgd_id.startswith('SGDID:'):
                    sgd_id = sgd_id[len('SGDID:'):]
                if not sgd_id:
                    sgd_id = next((
                        xref['gene_id'] for xref in xrefs
                        if xref['database'] == 'SGD'
                    ), '')
                ensembl_id = str(gene_info.get('ensembl_id', '') or '')
                if not ensembl_id:
                    ensembl_id = next((
                        xref['gene_id'] for xref in xrefs
                        if xref['database'] == 'Ensembl'
                    ), '')

                guide_entries = sorted({
                    item.strip()
                    for item in str(
                        gene_info.get('guide_targets', '') or ''
                    ).split(',')
                    if item.strip()
                })
                for guide_entry in guide_entries:
                    guide_identity = (pw_kfi_gene_key, guide_entry)
                    if guide_identity in seen_guides:
                        continue
                    seen_guides.add(guide_identity)
                    parts = guide_entry.split('/')
                    if len(parts) >= 6:
                        guide_info = {
                                            'pathway_id': pathway_id,
                                            'pathway_name': pathway_name,
                                            'pathway_gene_database': database,
                                            'pathway_gene_id': gene_id,
                                            'gene_key': gene_key,
                                            'pw_kfi_gene_key': pw_kfi_gene_key,
                                            'source_xref_keys': source_xref_keys,
                                            'source_databases': sorted({
                                                xref['database'] for xref in xrefs
                                            }),
                                            'annotation_source_database': gene_info.get(
                                                'annotation_source_database',
                                                'Unspecified',
                                            ),
                                            'matched_annotation_key': gene_info.get(
                                                'matched_annotation_key', ''
                                            ),
                                            'matched_annotation_id': gene_info.get(
                                                'matched_annotation_id', ''
                                            ),
                                            'resolved_annotation_subject_id': gene_info.get(
                                                'resolved_annotation_subject_id', ''
                                            ),
                                            'match_method': gene_info.get(
                                                'match_method', ''
                                            ),
                                            'sgd_id': sgd_id,
                                            'ensembl_id': ensembl_id,
                                            'unique_id': unique_id,
                                            'gene_name': gene_name,
                                            'description': description,
                                            'guide_unique_id': parts[0],
                                            'start_pos': parts[1],
                                            'end_pos': parts[2],
                                            'guide_sequence': parts[3],
                                            'strand': parts[7] if len(parts) > 7 else '',
                                            'MM0': parts[4],
                                            'MM1-2': parts[5],
                                            'off_target_match_score': parts[6] if len(parts) > 6 else ''
                        }
                        guides_data.append(guide_info)
            
            if guides_data:
                pathway_guide_data[pathway_id] = {
                    'pathway_name': pathway_name,
                    'guides': guides_data
                }
        
        self.logger.info(f"Pathway-specific guide CSV data generation complete.: {len(pathway_guide_data)} Pathway")
        return pathway_guide_data

    def generate_pathway_json_summary(self, pathways, pathway_guide_data=None, unique_id_to_sgd=None):
        """Generate a pathway summary JSON (including SVG data)."""
        # 統計情報計算
        total_pathways = len(pathways)
        fully_targetable = len([p for p in pathways if p['targetable_type'] == 'Fully targetable'])
        partially_targetable = len([p for p in pathways if p['targetable_type'] == 'Partially targetable'])
        not_targetable = len([p for p in pathways if p['targetable_type'] == 'Not targetable'])
        
        # パーセンテージ計算
        if total_pathways > 0:
            fully_targetable_percent = (fully_targetable / total_pathways) * 100
            partially_targetable_percent = (partially_targetable / total_pathways) * 100
            not_targetable_percent = (not_targetable / total_pathways) * 100
        else:
            fully_targetable_percent = partially_targetable_percent = not_targetable_percent = 0
        
        # チャートデータ生成
        chart_data = {
            'targetability': {
                'labels': ['Fully Targetable', 'Partially Targetable', 'Not Targetable'],
                'data': [fully_targetable, partially_targetable, not_targetable],
                'percentages': [fully_targetable_percent, partially_targetable_percent, not_targetable_percent],
                'colors': ['#28a745', '#ffc107', '#dc3545']
            }
        }
        
        # パスウェイSVGデータを事前生成
        pathway_svgs = {}
        svg_generation_stats = {
            'success': 0,
            'failed': 0,
            'no_gpml': 0,
            'failure_reasons': {},
        }
        
        self.logger.info(f"Started generating pathway SVG data.: {total_pathways} pathway")
        
        for i, pathway in enumerate(pathways):
            pathway_id = pathway['pathway_id']
            
            try:
                # 進捗表示
                if (i + 1) % 10 == 0 or i == 0:
                    self.logger.info(f"Pathway SVG generation progress: {i + 1}/{total_pathways}")
                
                # パスウェイ図データを生成
                diagram_data = self.get_pathway_diagram_data(pathway_id, pathways)

                if diagram_data:
                    pathway_svgs[pathway_id] = {
                        'pathway_name': diagram_data['pathway_name'],
                        'organism': diagram_data['organism'],
                        'svg_content': diagram_data['svg_content'],
                        'board_width': diagram_data['board_width'],
                        'board_height': diagram_data['board_height'],
                        'total_genes': diagram_data['total_genes'],
                        'evaluable_gene_nodes': diagram_data[
                            'evaluable_gene_nodes'
                        ],
                        'unique_xref_count': diagram_data[
                            'unique_xref_count'
                        ],
                        'total_gene_nodes': diagram_data['total_gene_nodes'],
                        'genes_with_guides': diagram_data['genes_with_guides'],
                        'gene_nodes_with_guides': diagram_data[
                            'gene_nodes_with_guides'
                        ],
                        'guide_status_map': diagram_data['guide_status_map'],
                        'unrenderable_gene_records': diagram_data[
                            'unrenderable_gene_records'
                        ],
                    }
                    svg_generation_stats['success'] += 1
                else:
                    # GPML file is missing or exact identifier data is unavailable.
                    gpml_file_path = pathway.get('gpml_file_path', '')
                    if not gpml_file_path or not os.path.exists(gpml_file_path):
                        svg_generation_stats['no_gpml'] += 1
                    else:
                        svg_generation_stats['failed'] += 1
                        svg_generation_stats['failure_reasons'][pathway_id] = (
                            'gpml_parse_or_exact_identifier_data_unavailable'
                        )

                    self.logger.debug(
                        f"Pathway {pathway_id} skipping SVG generation: "
                        "GPML or exact identifier issue"
                    )
                    
            except Exception as e:
                svg_generation_stats['failed'] += 1
                svg_generation_stats['failure_reasons'][pathway_id] = str(e)
                self.logger.warning(f"Pathway {pathway_id} SVG generation error: {e}")
        
        self.logger.info(f"Pathway SVG generation completed: Success={svg_generation_stats['success']}, "
                        f"Failure={svg_generation_stats['failed']}, No GPML={svg_generation_stats['no_gpml']}")
        
        summary_data = {
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'statistics': {
                'pathway_analysis': {
                    'total_pathways': total_pathways,
                    'fully_targetable_pathways': fully_targetable,
                    'partially_targetable_pathways': partially_targetable,
                    'not_targetable_pathways': not_targetable,
                    'fully_targetable_percent': fully_targetable_percent,
                    'partially_targetable_percent': partially_targetable_percent,
                    'not_targetable_percent': not_targetable_percent,
                    'chart_data': chart_data
                },
                'svg_generation': svg_generation_stats
            },
            'pathway_guide_data': pathway_guide_data if pathway_guide_data else {},
            'pathway_svgs': pathway_svgs
        }
        
        # JSON保存
        json_path = os.path.join(self.web_data_dir, 'pathway_summary.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(summary_data, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"Pathway summary JSON generation: {json_path}")
        return summary_data
    
    def create_pathway_html_template(self, pathways, summary_data, pathway_genes=None):
        """Generate a dedicated HTML template for the pathway."""
        # パスウェイデータを準備
        pathway_table_data = []
        for pathway in pathways:
            pathway_table_data.append({
                'pathway_id': pathway['pathway_id'],
                'pathway_name': pathway['pathway_name'],
                'total_genes': int(pathway.get('total_genes', 0)),
                'total_gene_nodes': int(
                    pathway.get('total_gene_nodes', 0) or 0
                ),
                'evaluable_gene_nodes': int(
                    pathway.get('evaluable_gene_node_count', 0) or 0
                ),
                'unique_xref_count': int(
                    pathway.get('unique_xref_count', 0) or 0
                ),
                'matched_genes': str(pathway.get('matched_genes', 0)),
                'genes_with_guides': int(pathway.get('genes_with_guides', 0)),
                'pw_kfi_score': float(pathway.get('pw_kfi_score', 0)),
                'targetable_type': pathway.get('targetable_type', 'Not targetable'),
                'coverage_rate': float(pathway.get('coverage_rate', 0))
            })
        
        # パスウェイ-遺伝子マッピングデータを準備
        pathway_genes_data = {}
        if pathway_genes:
            for pathway_id, genes in pathway_genes.items():
                pathway_genes_data[pathway_id] = []
                for gene in genes:
                    pathway_genes_data[pathway_id].append({
                        'pathway_gene_database': gene.get(
                            'pathway_gene_database', ''
                        ),
                        'pathway_gene_id': gene.get('pathway_gene_id', ''),
                        'gene_key': gene.get('gene_key', ''),
                        'pw_kfi_gene_key': gene.get('pw_kfi_gene_key', ''),
                        'pathway_xrefs': gene.get('pathway_xrefs', ''),
                        'source_databases': gene.get('source_databases', ''),
                        'node_occurrence_count': int(
                            gene.get('node_occurrence_count', 0) or 0
                        ),
                        'unique_xref_count': int(
                            gene.get('unique_xref_count', 0) or 0
                        ),
                        'canonical_subject_id': gene.get(
                            'canonical_subject_id', ''
                        ),
                        'annotation_source_database': gene.get(
                            'annotation_source_database', ''
                        ),
                        'matched_annotation_key': gene.get(
                            'matched_annotation_key', ''
                        ),
                        'matched_annotation_id': gene.get(
                            'matched_annotation_id', ''
                        ),
                        'resolved_annotation_subject_id': gene.get(
                            'resolved_annotation_subject_id', ''
                        ),
                        'match_method': gene.get('match_method', ''),
                        'sgd_id': gene.get('sgd_id', ''),
                        'ensembl_id': gene.get('ensembl_id', ''),
                        'gene_sgd_id': gene.get('gene_sgd_id', ''),
                        'unique_id': gene.get('unique_id', ''),
                        'gene_name': gene.get('gene_name', ''),
                        'has_guides': str(gene.get('has_guides', 'False')).lower() == 'true',
                        'guide_count': int(gene.get('guide_count', 0))
                    })
        
        # データをJSON形式で埋め込み
        embedded_data = {
            'summary': summary_data,
            'pathways': pathway_table_data,
            'pathway_genes': pathway_genes_data,
            'pathway_guide_data': summary_data.get('pathway_guide_data', {}),
            'pathway_svgs': summary_data.get('pathway_svgs', {})
        }
        
        html_template = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>tWhEEL - Pathway Analysis Report</title>
    <style>
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        
        html {{
            scroll-behavior: smooth;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.5;
            background-color: #f8f9fa;
            color: #212529;
        }}
        
        /* Layout */
        .report-layout {{
            display: flex;
            min-height: 100vh;
        }}
        
        /* Sidebar */
        .sidebar {{
            width: 280px;
            background-color: #fff;
            border-right: 1px solid #dee2e6;
            padding: 0;
            position: fixed;
            height: 100vh;
            overflow-y: auto;
            z-index: 1000;
        }}
        
        .sidebar-header {{
            padding: 1.5rem;
            border-bottom: 1px solid #dee2e6;
            text-align: center;
        }}
        
        .sidebar-header img {{
            max-width: 150px;
            height: auto;
            margin-bottom: 0.5rem;
        }}
        
        .sidebar-header h1 {{
            font-size: 1.75rem;
            font-weight: 600;
            color: #5cb85c;
            margin: 0;
        }}
        
        .sidebar-header .subtitle {{
            color: #6c757d;
            font-size: 0.85rem;
            margin: 0.5rem 0;
            font-weight: 400;
        }}
        
        .sidebar-header .meta {{
            color: #6c757d;
            font-size: 0.8rem;
            margin: 0.25rem 0 0 0;
        }}
        
        .sidebar-nav {{
            padding: 0;
        }}
        
        .sidebar-nav ul {{
            list-style: none;
            padding: 0;
            margin: 0;
        }}
        
        .sidebar-nav li {{
            margin: 0;
        }}
        
        .sidebar-nav a {{
            display: block;
            padding: 0.75rem 1.5rem;
            color: #495057;
            text-decoration: none;
            font-size: 0.9rem;
            border-left: 3px solid transparent;
            transition: all 0.2s ease;
        }}
        
        .sidebar-nav a:hover {{
            background-color: #f8f9fa;
            color: #5cb85c;
            border-left-color: #dee2e6;
        }}
        
        .sidebar-nav a.active {{
            background-color: #e8f5e9;
            color: #5cb85c;
            border-left-color: #5cb85c;
            font-weight: 500;
        }}
        
        /* Main content */
        .main-content {{
            margin-left: 280px;
            flex: 1;
            padding: 0;
        }}
        
        .content-header {{
            background-color: #fff;
            border-bottom: 1px solid #dee2e6;
            padding: 1.5rem 2rem;
        }}
        
        .content-header h1 {{
            font-size: 2rem;
            font-weight: 600;
            color: #212529;
            margin: 0;
        }}
        
        .content-header .subtitle {{
            color: #6c757d;
            font-size: 1rem;
            margin: 0.25rem 0 0 0;
        }}
        
        .content-body {{
            padding: 2rem;
        }}
        
        /* Sections */
        .content-section {{
            background-color: #fff;
            border: 1px solid #dee2e6;
            border-radius: 0.375rem;
            margin-bottom: 2rem;
            overflow: hidden;
        }}
        
        .section-header {{
            background-color: #f8f9fa;
            border-bottom: 1px solid #dee2e6;
            padding: 1rem 1.5rem;
        }}
        
        .section-header h2 {{
            font-size: 1.25rem;
            font-weight: 600;
            color: #212529;
            margin: 0;
        }}
        
        .section-body {{
            padding: 1.5rem;
        }}
        
        /* Stats Grid */
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 1.5rem;
            margin-bottom: 1.5rem;
        }}
        
        .stat-card {{
            background: white;
            padding: 1.5rem;
            border-radius: 8px;
            border: 1px solid #dee2e6;
            text-align: center;
        }}
        
        .stat-card.total {{ border-left: 4px solid #007bff; }}
        .stat-card.fully {{ border-left: 4px solid #28a745; }}
        .stat-card.partially {{ border-left: 4px solid #ffc107; }}
        .stat-card.not {{ border-left: 4px solid #dc3545; }}
        
        .stat-number {{
            font-size: 2rem;
            font-weight: bold;
            margin-bottom: 0.5rem;
        }}
        
        .stat-label {{
            font-size: 0.9rem;
            color: #666;
        }}
        
        .chart-container {{
            background: white;
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 1rem;
            border: 1px solid #dee2e6;
        }}
        
        .table-container {{
            background: white;
            overflow: hidden;
        }}
        
        .table-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 1rem 1.5rem;
            flex-wrap: wrap;
            gap: 1rem;
            border-bottom: 1px solid #dee2e6;
        }}
        
        .table-header h3 {{
            margin: 0;
            font-size: 1rem;
            font-weight: 600;
            color: #212529;
        }}
        
        .search-filter {{
            display: flex;
            gap: 1rem;
            align-items: center;
            flex-wrap: wrap;
        }}
        
        input[type="text"], select {{
            padding: 0.5rem;
            border: 1px solid #ced4da;
            border-radius: 4px;
            font-size: 0.9rem;
        }}
        
        .btn {{
            background: #5cb85c;
            color: white;
            border: none;
            padding: 0.5rem 1rem;
            border-radius: 4px;
            cursor: pointer;
            text-decoration: none;
            display: inline-flex;
            align-items: center;
            gap: 0.5rem;
            font-size: 0.9rem;
        }}
        
        .btn:hover {{
            background: #4cae4c;
        }}
        
        .btn.secondary {{
            background: #6c757d;
        }}
        
        .btn.secondary:hover {{
            background: #545b62;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
        }}
        
        th, td {{
            text-align: left;
            padding: 0.75rem 1rem;
            border-bottom: 1px solid #dee2e6;
        }}
        
        th {{
            background-color: #f8f9fa;
            font-weight: 600;
            color: #495057;
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        
        tr:hover {{
            background-color: #f8f9fa;
        }}
        
        .targetable-fully {{ color: #28a745; font-weight: 600; }}
        .targetable-partially {{ color: #856404; font-weight: 600; }}
        .targetable-not {{ color: #dc3545; font-weight: 600; }}
        
        .modal {{
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.5);
        }}
        
        .modal-content {{
            background-color: white;
            margin: 5% auto;
            padding: 2rem;
            border-radius: 8px;
            width: 90%;
            max-width: 1000px;
            max-height: 80vh;
            overflow-y: auto;
        }}
        
        .close {{
            color: #aaa;
            float: right;
            font-size: 28px;
            font-weight: bold;
            cursor: pointer;
        }}
        
        .close:hover {{
            color: black;
        }}
        
        .pagination {{
            display: flex;
            justify-content: center;
            align-items: center;
            padding: 1rem;
            gap: 0.5rem;
        }}
        
        .page-btn {{
            padding: 0.5rem 0.75rem;
            border: 1px solid #dee2e6;
            background: white;
            cursor: pointer;
            border-radius: 4px;
        }}
        
        .page-btn:hover, .page-btn.active {{
            background: #5cb85c;
            color: white;
            border-color: #5cb85c;
        }}
        
        .loading {{
            text-align: center;
            padding: 2rem;
            color: #666;
        }}
        
        @media (max-width: 768px) {{
            .sidebar {{
                display: none;
            }}
            
            .main-content {{
                margin-left: 0;
            }}
            
            .stats-grid {{
                grid-template-columns: 1fr;
            }}
            
            .search-filter {{
                flex-direction: column;
                align-items: stretch;
            }}
            
            table {{
                font-size: 0.8rem;
            }}
            
            th, td {{
                padding: 0.5rem;
            }}
        }}
    </style>
</head>
<body>
    <!-- データ埋め込み -->
    <script type="application/json" id="pathway-data">{json.dumps(embedded_data, ensure_ascii=False, indent=2)}</script>
    
    <div class="report-layout">
        <!-- Sidebar -->
        <aside class="sidebar">
            <div class="sidebar-header">
                <img src="static/img/logo.png" alt="tWhEEL Logo">
                <p class="subtitle">A pipeline to aggregate results from gene functional analysis, guide RNA design, and pathway-based knock-out feasibility assessment into a single report.</p>
                <p class="meta">Report generated on {summary_data.get('timestamp', '')}</p>
            </div>
            <nav class="sidebar-nav">
                <ul>
                    <li><a href="#summary-stats" class="nav-link active">Summary Statistics</a></li>
                    <li><a href="#pathway-analysis" class="nav-link">Pathway Analysis</a></li>
                </ul>
            </nav>
        </aside>
        
        <!-- Main Content -->
        <main class="main-content">
            <header class="content-header">
                <h1>tWhEEL</h1>
                <p class="subtitle">Transcriptome Workflow for Hierarchical Exploration and Expression Levels</p>
            </header>
            <div class="content-body">
                <section id="summary-stats" class="content-section">
                    <div class="section-header">
                        <h2>Summary Statistics</h2>
                    </div>
                    <div class="section-body">
        <div class="stats-grid">
            <div class="stat-card total">
                <div class="stat-number" id="total-pathways">-</div>
                <div class="stat-label">Total Pathways</div>
            </div>
            <div class="stat-card fully">
                <div class="stat-number" id="fully-targetable">-</div>
                <div class="stat-label">Fully Targetable</div>
            </div>
            <div class="stat-card partially">
                <div class="stat-number" id="partially-targetable">-</div>
                <div class="stat-label">Partially Targetable</div>
            </div>
            <div class="stat-card not">
                <div class="stat-number" id="not-targetable">-</div>
                <div class="stat-label">Not Targetable</div>
            </div>
        </div>
        
        <div class="chart-container">
            <h3>Pathway Targetability Distribution</h3>
            <div id="targetabilityChart" style="display: flex; justify-content: center; align-items: center;"></div>
        </div>
                    </div>
                </section>
                
                <section id="pathway-analysis" class="content-section">
                    <div class="section-header">
                        <h2>Pathway Analysis</h2>
                    </div>
                    <div class="section-body" style="padding: 0;">
        <div class="table-container">
            <div class="table-header">
                <h3>Pathway Analysis Results</h3>
                <div class="search-filter">
                    <input type="text" id="search-input" placeholder="Search pathways...">
                    <select id="targetable-filter">
                        <option value="">All Types</option>
                        <option value="Fully targetable">Fully Targetable</option>
                        <option value="Partially targetable">Partially Targetable</option>
                        <option value="Not targetable">Not Targetable</option>
                    </select>
                    <button class="btn" onclick="downloadPathwayCSV()">
                        📥 Download CSV
                    </button>
                </div>
            </div>
            
            <div id="table-loading" class="loading">Loading pathway data...</div>
            
            <table id="pathway-table" style="display: none;">
                <thead>
                    <tr>
                        <th onclick="sortTable(0)">Pathway ID</th>
                        <th onclick="sortTable(1)">Pathway Name</th>
                        <th onclick="sortTable(2)">Raw GPML Nodes</th>
                        <th onclick="sortTable(3)">Evaluable Nodes</th>
                        <th onclick="sortTable(4)">Unique Xrefs</th>
                        <th onclick="sortTable(5)">Unique PW-KFI Genes</th>
                        <th onclick="sortTable(6)">PW-KFI Genes with Guides</th>
                        <th onclick="sortTable(7)">PW-KFI Score</th>
                        <th onclick="sortTable(8)">Targetability</th>
                        <th onclick="sortTable(9)">Annotation Match Coverage</th>
                        <th>Actions</th>
                    </tr>
                </thead>
                <tbody id="pathway-table-body">
                </tbody>
            </table>
            
            <div class="pagination" id="pagination" style="display: none;">
            </div>
        </div>
                    </div>
                </section>
            </div>
        </main>
    </div>
    
    <!-- パスウェイ遺伝子詳細モーダル -->
    <div id="pathway-genes-modal" class="modal">
        <div class="modal-content" style="max-width: 90%; width: 90%;">
            <span class="close" onclick="closeModal('pathway-genes-modal')">&times;</span>
            <h3 id="modal-pathway-name">Pathway Genes</h3>
            <div class="table-container" style="max-height: 70vh; overflow-y: auto;">
                <table id="pathway-genes-table">
                    <thead>
                        <tr>
                            <th>Gene Name</th>
                            <th>Pathway Xref Database</th>
                            <th>Pathway Xref ID</th>
                            <th>Node Occurrences</th>
                            <th>Annotation Source</th>
                            <th>Description</th>
                            <th>Match Status</th>
                            <th>Failure Reason</th>
                            <th>Guide Status</th>
                            <th>Verified Link</th>
                        </tr>
                    </thead>
                    <tbody id="pathway-genes-table-body">
                    </tbody>
                </table>
                <div id="pathway-genes-empty" style="display: none; text-align: center; padding: 2rem; color: #666;">
                    This pathway does not contain genetic information.
                </div>
            </div>
        </div>
    </div>
    
    <!-- パスウェイ図表示モーダル -->
    <div id="pathway-diagram-modal" class="modal">
        <div class="modal-content" style="max-width: 95%; width: 95%;">
            <span class="close" onclick="closeModal('pathway-diagram-modal')">&times;</span>
            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 1rem;">
                <h3 id="diagram-pathway-name">Pathway Diagram</h3>
                <div>
                    <button class="btn secondary" onclick="zoomDiagram(0.8)">🔍-</button>
                    <button class="btn secondary" onclick="zoomDiagram(1.25)">🔍+</button>
                    <button class="btn" onclick="resetZoom()">Reset</button>
                </div>
            </div>
            <div id="pathway-diagram-loading" class="loading">Loading pathway diagram......</div>
            <div id="pathway-diagram-container" style="display: none; border: 1px solid #ddd; overflow: auto; max-height: 70vh; background: white;">
                <div id="pathway-diagram-content"></div>
            </div>
            <div id="pathway-diagram-info" style="margin-top: 1rem; font-size: 0.9rem; color: #666;"></div>
        </div>
    </div>
    
    <!-- 遺伝子詳細情報モーダル -->
    <div id="gene-detail-modal" class="modal">
        <div class="modal-content" style="max-width: 90%; width: 90%;">
            <span class="close" onclick="closeModal('gene-detail-modal')">&times;</span>
            <h3 id="gene-detail-title">Gene Details</h3>
            <div id="gene-basic-info" style="background: #f8f9fa; padding: 1rem; margin-bottom: 1rem; border-radius: 4px;">
                <!-- Basic genetic information is displayed here. -->
            </div>
            <div class="table-container" style="max-height: 60vh; overflow-y: auto;">
                <h4>Guide RNA Details</h4>
                <table id="gene-guides-table">
                    <thead>
                        <tr>
                            <th>Guide Sequence</th>
                            <th>Start Position</th>
                            <th>End Position</th>
                            <th>MM0</th>
                            <th>MM1-2</th>
                            <th>Off-target Match Score</th>
                            <th>Strand</th>
                        </tr>
                    </thead>
                    <tbody id="gene-guides-table-body">
                    </tbody>
                </table>
                <div id="gene-guides-empty" style="display: none; text-align: center; padding: 2rem; color: #666;">
                    This gene does not contain guide RNA design information.
                </div>
            </div>
        </div>
    </div>
    
    <script>
        let pathwayData = [];
        let filteredData = [];
        let currentPage = 1;
        const itemsPerPage = 25;
        let sortColumn = 0;
        let sortDirection = 'asc';
        
        // 初期化
        document.addEventListener('DOMContentLoaded', function() {{
            loadPathwayData();
            setupEventListeners();
        }});
        
        function setupEventListeners() {{
            document.getElementById('search-input').addEventListener('input', filterTable);
            document.getElementById('targetable-filter').addEventListener('change', filterTable);
        }}
        
        function loadPathwayData() {{
            try {{
                // 埋め込まれたデータを読み込み
                const dataElement = document.getElementById('pathway-data');
                if (!dataElement) {{
                    throw new Error('Embedded data not found.');
                }}
                
                const data = JSON.parse(dataElement.textContent);
                
                // サマリーデータを表示
                loadPathwaySummary(data.summary);
                
                // パスウェイテーブルデータをロード
                pathwayData = data.pathways;
                filteredData = [...pathwayData];
                
                document.getElementById('table-loading').style.display = 'none';
                document.getElementById('pathway-table').style.display = 'table';
                document.getElementById('pagination').style.display = 'flex';
                
                renderTable();
                
            }} catch (error) {{
                console.error('Error loading data:', error);
                document.getElementById('table-loading').textContent = 'Failed to load data.: ' + error.message;
            }}
        }}
        
        function loadPathwaySummary(summaryData) {{
            try {{
                const stats = summaryData.statistics.pathway_analysis;
                document.getElementById('total-pathways').textContent = stats.total_pathways;
                document.getElementById('fully-targetable').textContent = stats.fully_targetable_pathways;
                document.getElementById('partially-targetable').textContent = stats.partially_targetable_pathways;
                document.getElementById('not-targetable').textContent = stats.not_targetable_pathways;
                
                // チャート生成
                createTargetabilityChart(stats.chart_data.targetability);
                
            }} catch (error) {{
                console.error('Summary data loading error:', error);
            }}
        }}
        
        function createTargetabilityChart(chartData) {{
            const chartContainer = document.getElementById('targetabilityChart');
            if (!chartContainer) return;
            
            // Create horizontal bar chart using SVG
            const chartSVG = createPathwayTargetabilityBarChart(chartData);
            
            chartContainer.innerHTML = '';
            chartContainer.appendChild(chartSVG);
            
            console.log('Pathway Targetability bar chart rendered');
        }}
        
        function createPathwayTargetabilityBarChart(chartData) {{
            const width = 900;
            const height = 300;
            const margin = {{ top: 80, right: 50, bottom: 70, left: 80 }};
            
            // Create SVG element
            const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
            svg.setAttribute('width', width);
            svg.setAttribute('height', height);
            svg.style.fontFamily = 'Arial, sans-serif';
            svg.style.fontSize = '12px';
            
            // Chart dimensions
            const chartWidth = width - margin.left - margin.right;
            const chartHeight = height - margin.top - margin.bottom;
            
            // Bar position and height
            const barY = margin.top + (chartHeight - 60) / 2;
            const barHeight = 60;
            
            // Calculate segment widths based on percentages
            const fullyPercent = chartData.percentages[0] || 0;
            const partiallyPercent = chartData.percentages[1] || 0;
            const notPercent = chartData.percentages[2] || 0;
            
            // Total count for statistics
            const totalPathways = chartData.data.reduce((sum, count) => sum + count, 0);
            
            // Segment data with proper positioning
            const segments = [
                {{
                    label: 'Fully Targetable',
                    percent: fullyPercent,
                    count: chartData.data[0],
                    color: '#28a745',
                    x: 0
                }},
                {{
                    label: 'Partially Targetable',
                    percent: partiallyPercent,
                    count: chartData.data[1],
                    color: '#ffc107',
                    x: (fullyPercent / 100) * chartWidth
                }},
                {{
                    label: 'Not Targetable',
                    percent: notPercent,
                    count: chartData.data[2],
                    color: '#dc3545',
                    x: ((fullyPercent + partiallyPercent) / 100) * chartWidth
                }}
            ];
            
            // Title
            const title = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            title.setAttribute('x', width / 2);
            title.setAttribute('y', 20);
            title.setAttribute('text-anchor', 'middle');
            title.setAttribute('fill', '#333');
            title.setAttribute('font-weight', 'bold');
            title.setAttribute('font-size', '16px');
            title.textContent = 'Pathway Targetability Distribution';
            svg.appendChild(title);
            
            // Legend - arranged in single row
            const legendStartX = 150;
            const legendY = 45;
            const legendSquareSize = 12;
            const legendSpacingX = 200;
            
            segments.forEach((segment, i) => {{
                const legendX = legendStartX + (i * legendSpacingX);
                
                // Legend square
                const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                rect.setAttribute('x', legendX);
                rect.setAttribute('y', legendY - legendSquareSize / 2);
                rect.setAttribute('width', legendSquareSize);
                rect.setAttribute('height', legendSquareSize);
                rect.setAttribute('fill', segment.color);
                svg.appendChild(rect);
                
                // Legend text
                const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                text.setAttribute('x', legendX + legendSquareSize + 5);
                text.setAttribute('y', legendY + 3);
                text.setAttribute('fill', '#333');
                text.setAttribute('font-size', '11px');
                text.textContent = segment.label;
                svg.appendChild(text);
            }});
            
            // Y-axis label
            const yLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            yLabel.setAttribute('x', 25);
            yLabel.setAttribute('y', barY + barHeight / 2);
            yLabel.setAttribute('text-anchor', 'middle');
            yLabel.setAttribute('fill', '#333');
            yLabel.setAttribute('font-weight', 'bold');
            yLabel.textContent = 'Pathways';
            svg.appendChild(yLabel);
            
            // Draw stacked bar segments
            segments.forEach(segment => {{
                if (segment.percent > 0) {{
                    const segmentWidth = (segment.percent / 100) * chartWidth;
                    
                    // Segment rectangle
                    const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
                    rect.setAttribute('x', margin.left + segment.x);
                    rect.setAttribute('y', barY);
                    rect.setAttribute('width', segmentWidth);
                    rect.setAttribute('height', barHeight);
                    rect.setAttribute('fill', segment.color);
                    rect.setAttribute('stroke', '#fff');
                    rect.setAttribute('stroke-width', '1');
                    svg.appendChild(rect);
                }}
            }});
            
            // Combined statistics text
            const statsText = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            statsText.setAttribute('x', margin.left + chartWidth / 2);
            statsText.setAttribute('y', barY + barHeight / 2);
            statsText.setAttribute('text-anchor', 'middle');
            statsText.setAttribute('fill', '#333');
            statsText.setAttribute('font-weight', 'bold');
            statsText.setAttribute('font-size', '12px');
            statsText.textContent = `F:${{segments[0].count}} P:${{segments[1].count}} N:${{segments[2].count}} n:${{totalPathways}}`;
            svg.appendChild(statsText);
            
            // X-axis scale
            const scaleY = barY + barHeight + 20;
            const scaleMarks = [0, 20, 40, 60, 80, 100];
            
            // X-axis line
            const xAxisLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
            xAxisLine.setAttribute('x1', margin.left);
            xAxisLine.setAttribute('y1', scaleY);
            xAxisLine.setAttribute('x2', margin.left + chartWidth);
            xAxisLine.setAttribute('y2', scaleY);
            xAxisLine.setAttribute('stroke', '#333');
            xAxisLine.setAttribute('stroke-width', '1');
            svg.appendChild(xAxisLine);
            
            // X-axis scale marks and labels
            scaleMarks.forEach(mark => {{
                const x = margin.left + (mark / 100) * chartWidth;
                
                // Scale mark
                const tick = document.createElementNS('http://www.w3.org/2000/svg', 'line');
                tick.setAttribute('x1', x);
                tick.setAttribute('y1', scaleY);
                tick.setAttribute('x2', x);
                tick.setAttribute('y2', scaleY + 5);
                tick.setAttribute('stroke', '#333');
                tick.setAttribute('stroke-width', '1');
                svg.appendChild(tick);
                
                // Scale label
                const label = document.createElementNS('http://www.w3.org/2000/svg', 'text');
                label.setAttribute('x', x);
                label.setAttribute('y', scaleY + 18);
                label.setAttribute('text-anchor', 'middle');
                label.setAttribute('fill', '#333');
                label.setAttribute('font-size', '10px');
                label.textContent = mark;
                svg.appendChild(label);
            }});
            
            // X-axis label
            const xLabel = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            xLabel.setAttribute('x', margin.left + chartWidth / 2);
            xLabel.setAttribute('y', height - 10);
            xLabel.setAttribute('text-anchor', 'middle');
            xLabel.setAttribute('fill', '#333');
            xLabel.setAttribute('font-weight', 'bold');
            xLabel.textContent = '% PW-KFI';
            svg.appendChild(xLabel);
            
            return svg;
        }}
        
        function renderTable() {{
            const tbody = document.getElementById('pathway-table-body');
            const startIndex = (currentPage - 1) * itemsPerPage;
            const endIndex = startIndex + itemsPerPage;
            const paginatedData = filteredData.slice(startIndex, endIndex);
            
            tbody.innerHTML = '';
            
            paginatedData.forEach(pathway => {{
                const row = tbody.insertRow();
                const targetableClass = pathway.targetable_type.toLowerCase().replace(' ', '-');
                
                row.innerHTML = `
                    <td>${{pathway.pathway_id}}</td>
                    <td>${{pathway.pathway_name}}</td>
                    <td>${{pathway.total_gene_nodes}}</td>
                    <td>${{pathway.evaluable_gene_nodes}}</td>
                    <td>${{pathway.unique_xref_count}}</td>
                    <td>${{pathway.total_genes}}</td>
                    <td>${{pathway.genes_with_guides}}</td>
                    <td>${{pathway.pw_kfi_score.toFixed(3)}}</td>
                    <td class="targetable-${{targetableClass}}">${{pathway.targetable_type}}</td>
                    <td>${{(pathway.coverage_rate * 100).toFixed(1)}}%</td>
                    <td>
                        <button class="btn secondary" onclick="showPathwayGenes('${{pathway.pathway_id}}')">View Genes</button>
                        <button class="btn" onclick="showPathwayDiagram('${{pathway.pathway_id}}')">🔍 View Pathway</button>
                        <button class="btn" onclick="downloadPathwayGuides('${{pathway.pathway_id}}')">📥 Guides</button>
                    </td>
                `;
            }});
            
            renderPagination();
        }}
        
        function renderPagination() {{
            const totalPages = Math.ceil(filteredData.length / itemsPerPage);
            const pagination = document.getElementById('pagination');
            pagination.innerHTML = '';
            
            // Previous button
            if (currentPage > 1) {{
                const prevBtn = document.createElement('button');
                prevBtn.className = 'page-btn';
                prevBtn.textContent = 'Previous';
                prevBtn.onclick = () => {{ currentPage--; renderTable(); }};
                pagination.appendChild(prevBtn);
            }}
            
            // Page numbers
            for (let i = 1; i <= totalPages; i++) {{
                if (i === currentPage || i === 1 || i === totalPages || Math.abs(i - currentPage) <= 2) {{
                    const pageBtn = document.createElement('button');
                    pageBtn.className = 'page-btn' + (i === currentPage ? ' active' : '');
                    pageBtn.textContent = i;
                    pageBtn.onclick = () => {{ currentPage = i; renderTable(); }};
                    pagination.appendChild(pageBtn);
                }}
            }}
            
            // Next button
            if (currentPage < totalPages) {{
                const nextBtn = document.createElement('button');
                nextBtn.className = 'page-btn';
                nextBtn.textContent = 'Next';
                nextBtn.onclick = () => {{ currentPage++; renderTable(); }};
                pagination.appendChild(nextBtn);
            }}
        }}
        
        function filterTable() {{
            const searchTerm = document.getElementById('search-input').value.toLowerCase();
            const targetableFilter = document.getElementById('targetable-filter').value;
            
            filteredData = pathwayData.filter(pathway => {{
                const matchesSearch = !searchTerm || 
                    pathway.pathway_id.toLowerCase().includes(searchTerm) ||
                    pathway.pathway_name.toLowerCase().includes(searchTerm);
                
                const matchesTargetable = !targetableFilter || pathway.targetable_type === targetableFilter;
                
                return matchesSearch && matchesTargetable;
            }});
            
            currentPage = 1;
            renderTable();
        }}
        
        function sortTable(column) {{
            if (sortColumn === column) {{
                sortDirection = sortDirection === 'asc' ? 'desc' : 'asc';
            }} else {{
                sortColumn = column;
                sortDirection = 'asc';
            }}
            
            const columns = [
                'pathway_id', 'pathway_name', 'total_gene_nodes',
                'evaluable_gene_nodes', 'unique_xref_count', 'total_genes',
                'genes_with_guides', 'pw_kfi_score', 'targetable_type',
                'coverage_rate'
            ];
            const sortKey = columns[column];
            
            filteredData.sort((a, b) => {{
                let valueA = a[sortKey];
                let valueB = b[sortKey];
                
                if (typeof valueA === 'string') {{
                    valueA = valueA.toLowerCase();
                    valueB = valueB.toLowerCase();
                }}
                
                if (sortDirection === 'asc') {{
                    return valueA < valueB ? -1 : valueA > valueB ? 1 : 0;
                }} else {{
                    return valueA > valueB ? -1 : valueA < valueB ? 1 : 0;
                }}
            }});
            
            renderTable();
        }}
        
        function showPathwayGenes(pathwayId) {{
            try {{
                // Retrieve pathway SVG information from embedded data
                const dataElement = document.getElementById('pathway-data');
                const data = JSON.parse(dataElement.textContent);
                const pathwaySvgData = data.pathway_svgs[pathwayId];
                
                if (!pathwaySvgData || !pathwaySvgData.guide_status_map) {{
                    alert('Genetic information for this pathway could not be found.: ' + pathwayId);
                    return;
                }}
                
                // パスウェイ名を取得
                const pathway = data.pathways.find(p => p.pathway_id === pathwayId);
                const pathwayName = pathway ? pathway.pathway_name : pathwayId;
                
                // モーダルタイトルを設定
                document.getElementById('modal-pathway-name').textContent = pathwayName + ' (' + pathwayId + ') - Genes';
                
                // 遺伝子基本情報テーブルを構築
                const tbody = document.getElementById('pathway-genes-table-body');
                const emptyDiv = document.getElementById('pathway-genes-empty');
                tbody.innerHTML = '';
                
                // guide_status_mapから遺伝子情報を取得
                const guideStatusMap = pathwaySvgData.guide_status_map;
                const geneList = Object.keys(guideStatusMap);
                
                if (geneList.length === 0) {{
                    // 遺伝子情報がない場合
                    emptyDiv.style.display = 'block';
                    tbody.parentElement.parentElement.style.display = 'none';
                }} else {{
                    // 遺伝子基本情報を表示
                    emptyDiv.style.display = 'none';
                    tbody.parentElement.parentElement.style.display = 'block';
                    
                    geneList.forEach(geneKey => {{
                        const geneInfo = guideStatusMap[geneKey];
                        const hasGuides = geneInfo.has_guides;
                        const geneName = geneInfo.gene_name || 'Unknown';
                        const pathwayDatabase = geneInfo.pathway_gene_database || 'Unavailable';
                        const pathwayGeneId = geneInfo.pathway_gene_id || 'Unavailable';
                        const nodeOccurrenceCount = geneInfo.node_occurrence_count || 0;
                        const annotationSource = geneInfo.annotation_source_database || 'Unspecified';
                        const matchStatus = geneInfo.match_status || 'unavailable';
                        const failureReason = geneInfo.failure_reason || '';
                        
                        // Descriptionを取得（pathway_guide_dataから）
                        let description = geneInfo.description || 'N/A';
                        const pathwayGuideData = data.pathway_guide_data[pathwayId];
                        if (pathwayGuideData && pathwayGuideData.guides) {{
                            const guideInfo = pathwayGuideData.guides.find(g => g.gene_key === geneKey);
                            if (guideInfo && guideInfo.description) {{
                                description = guideInfo.description;
                            }}
                        }}
                        
                        const row = tbody.insertRow();
                        row.innerHTML = `
                            <td><strong>${{geneName}}</strong></td>
                            <td>${{pathwayDatabase}}</td>
                            <td style="font-family: monospace; background: #f8f9fa;">${{pathwayGeneId}}</td>
                            <td>${{nodeOccurrenceCount}}</td>
                            <td>${{annotationSource}}</td>
                            <td title="${{description}}">${{
                                description && description.length > 40 ? 
                                description.substring(0, 40) + '...' : 
                                description
                            }}</td>
                            <td>${{matchStatus}}</td>
                            <td style="font-family: monospace;">${{failureReason || '—'}}</td>
                            <td>
                                ${{matchStatus === 'matched' ? `
                                    <span style="color: ${{hasGuides ? '#28a745' : '#dc3545'}}; font-weight: bold;">
                                        ${{hasGuides ? '✓ Available' : '✗ Not Available'}}
                                    </span>
                                    ${{hasGuides ? ` (${{geneInfo.guide_count}} guides)` : ''}}
                                ` : '<span style="color: #6c757d; font-weight: bold;">Unavailable</span>'}}
                            </td>
                            <td>
                                ${{geneInfo.sgd_id ? `
                                    <button class="btn secondary" onclick="openSgdLink('${{geneInfo.sgd_id}}')"
                                        title="Open verified SGD annotation">
                                        🔗 View SGD annotation
                                    </button>
                                ` : '<span title="No verified external URL is configured">Unavailable</span>'}}
                            </td>
                        `;
                    }});
                    
                    // 統計情報を追加
                    const statsRow = tbody.insertRow();
                    statsRow.style.backgroundColor = '#e9ecef';
                    statsRow.style.fontWeight = 'bold';
                    const geneNodesWithGuides = geneList.reduce((total, geneKey) =>
                        total + (guideStatusMap[geneKey].has_guides
                            ? (guideStatusMap[geneKey].node_occurrence_count || 0)
                            : 0), 0
                    );
                    const totalGeneNodes = pathwaySvgData.total_gene_nodes || 0;
                    const evaluableGeneNodes = pathwaySvgData.evaluable_gene_nodes || 0;
                    const unrenderableRecords = pathwaySvgData.unrenderable_gene_records || [];
                    const unrenderableReasons = [...new Set(
                        unrenderableRecords.map(record => record.failure_reason || 'unavailable')
                    )].join(', ');
                    const totalGenes = pathwaySvgData.total_genes || 0;
                    const genesWithGuides = pathwaySvgData.genes_with_guides || 0;
                    const uniqueXrefCount = pathwaySvgData.unique_xref_count || 0;
                    const coverage = totalGenes > 0
                        ? ((genesWithGuides / totalGenes) * 100).toFixed(1)
                        : '0.0';
                    
                    statsRow.innerHTML = `
                        <td colspan="10" style="text-align: center; padding: 1rem;">
                            <strong>統計情報:</strong> 
                            ${{totalGeneNodes}} total GPML gene nodes,
                            ${{evaluableGeneNodes}} evaluable nodes,
                            ${{uniqueXrefCount}} unique Database::ID xrefs,
                            ${{totalGenes}} unique PW-KFI genes,
                            ${{genesWithGuides}} PW-KFI genes with guide RNA targets,
                            ${{geneNodesWithGuides}} colored GPML nodes
                            (${{coverage}}% PW-KFI gene coverage)
                            ${{unrenderableRecords.length > 0
                                ? `<br>${{unrenderableRecords.length}} nodes cannot be keyed for diagram lookup: ${{unrenderableReasons}}`
                                : ''}}
                        </td>
                    `;
                }}
                
                // モーダルを表示
                document.getElementById('pathway-genes-modal').style.display = 'block';
                
            }} catch (error) {{
                console.error('Genetic data display error:', error);
                alert('Failed to display genetic data.: ' + error.message);
            }}
        }}
        
        function openSgdLink(sgdId) {{
            if (sgdId && sgdId !== 'N/A') {{
                const sgdUrl = `https://www.yeastgenome.org/locus/${{sgdId}}`;
                window.open(sgdUrl, '_blank');
            }} else {{
                alert('No valid SGD ID was found.');
            }}
        }}
        
        function copyToClipboard(text) {{
            // クリップボードにテキストをコピー
            if (navigator.clipboard && window.isSecureContext) {{
                navigator.clipboard.writeText(text).then(() => {{
                    // 成功時のフィードバック（簡易的）
                    const originalTitle = event.target.title;
                    event.target.title = 'Copied!';
                    setTimeout(() => {{
                        event.target.title = originalTitle;
                    }}, 1000);
                }}).catch(err => {{
                    console.error('Clipboard copy error:', err);
                }});
            }} else {{
                // フォールバック: 古いブラウザ対応
                const textArea = document.createElement('textarea');
                textArea.value = text;
                textArea.style.position = 'fixed';
                textArea.style.opacity = '0';
                document.body.appendChild(textArea);
                textArea.focus();
                textArea.select();
                try {{
                    document.execCommand('copy');
                    const originalTitle = event.target.title;
                    event.target.title = 'Copied!';
                    setTimeout(() => {{
                        event.target.title = originalTitle;
                    }}, 1000);
                }} catch (err) {{
                    console.error('Clipboard copy error:', err);
                }}
                document.body.removeChild(textArea);
            }}
        }}
        
        function closeModal(modalId) {{
            if (modalId) {{
                document.getElementById(modalId).style.display = 'none';
            }} else {{
                // 後方互換性のため
                document.getElementById('pathway-genes-modal').style.display = 'none';
            }}
        }}
        
        // パスウェイ図表示機能の変数
        let currentZoom = 1.0;
        let currentPathwayId = null;
        
        function showPathwayDiagram(pathwayId) {{
            try {{
                currentPathwayId = pathwayId;
                
                // パスウェイ名を取得
                const pathway = pathwayData.find(p => p.pathway_id === pathwayId);
                if (!pathway) {{
                    alert('Pathway not found.: ' + pathwayId);
                    return;
                }}
                
                // モーダルタイトルを設定
                document.getElementById('diagram-pathway-name').textContent = 
                    pathway.pathway_name + ' (' + pathwayId + ')';
                
                // ローディング表示
                document.getElementById('pathway-diagram-loading').style.display = 'block';
                document.getElementById('pathway-diagram-container').style.display = 'none';
                document.getElementById('pathway-diagram-info').textContent = '';
                
                // モーダルを表示
                document.getElementById('pathway-diagram-modal').style.display = 'block';
                
                // パスウェイ図を非同期で読み込み
                loadPathwayDiagram(pathwayId, pathway);
                
            }} catch (error) {{
                console.error('Pathway diagram display error:', error);
                alert('Failed to display the pathway diagram.: ' + error.message);
            }}
        }}
        
        function loadPathwayDiagram(pathwayId, pathway) {{
            // 埋め込みデータから実際のパスウェイSVGデータを取得
            setTimeout(() => {{
                try {{
                    const dataElement = document.getElementById('pathway-data');
                    const data = JSON.parse(dataElement.textContent);
                    const pathwaySvgData = data.pathway_svgs[pathwayId];
                    
                    if (!pathwaySvgData) {{
                        // SVGデータが見つからない場合のエラー表示
                        document.getElementById('pathway-diagram-loading').innerHTML = 
                            `<div style="color: red;">The diagram data for this pathway (${{pathwayId}}) could not be found.<br>` +
                            `The GPML file may not exist, or the parsing may have failed.</div>`;
                        return;
                    }}
                    
                    // 実際のSVGデータを表示
                    document.getElementById('pathway-diagram-content').innerHTML = pathwaySvgData.svg_content;
                    document.getElementById('pathway-diagram-loading').style.display = 'none';
                    document.getElementById('pathway-diagram-container').style.display = 'block';
                    
                    // 実際の統計情報を表示
                    const totalGeneNodes = pathwaySvgData.total_gene_nodes || 0;
                    const evaluableGeneNodes = pathwaySvgData.evaluable_gene_nodes || 0;
                    const uniqueXrefCount = pathwaySvgData.unique_xref_count || 0;
                    const totalGenes = pathwaySvgData.total_genes || 0;
                    const genesWithGuides = pathwaySvgData.genes_with_guides || 0;
                    const geneNodesWithGuides = pathwaySvgData.gene_nodes_with_guides || 0;
                    const genesWithoutGuides = totalGenes - genesWithGuides;
                    const coverage = totalGenes > 0
                        ? ((genesWithGuides / totalGenes) * 100).toFixed(1)
                        : '0.0';
                    
                    document.getElementById('pathway-diagram-info').innerHTML = `
                        <strong>Pathway:</strong> ${{pathwaySvgData.pathway_name}} | 
                        <strong>Total GPML gene nodes:</strong> ${{totalGeneNodes}} | 
                        <strong>Evaluable nodes:</strong> ${{evaluableGeneNodes}} | 
                        <strong>Unique Database::ID xrefs:</strong> ${{uniqueXrefCount}} | 
                        <strong>Unique PW-KFI genes:</strong> ${{totalGenes}} | 
                        <strong>PW-KFI genes with guides:</strong> ${{genesWithGuides}} | 
                        <strong>PW-KFI genes without guides:</strong> ${{genesWithoutGuides}} | 
                        <strong>Colored GPML nodes:</strong> ${{geneNodesWithGuides}} | 
                        <strong>PW-KFI gene coverage:</strong> ${{coverage}}%
                    `;
                    
                    // SVG要素にクリックイベントを追加
                    addGeneClickEvents();
                    
                    // ズームをリセット
                    currentZoom = 1.0;
                    updateDiagramZoom();
                    
                }} catch (error) {{
                    console.error('パスウェイ図表示エラー:', error);
                    document.getElementById('pathway-diagram-loading').innerHTML = 
                        '<div style="color: red;">Failed to display the pathway diagram.: ' + error.message + '</div>';
                }}
            }}, 500); // 短い遅延で実際のデータを表示
        }}
        
        function generateSampleSVG(pathwayId, genes) {{
            const width = 600;
            const height = 500;
            
            let svg = `<svg width="${{width}}" height="${{height}}" viewBox="0 0 ${{width}} ${{height}}" xmlns="http://www.w3.org/2000/svg">`;
            svg += '<defs><style>';
            svg += '.gene-node {{ stroke-width: 2; cursor: pointer; }}';
            svg += '.gene-with-guides {{ fill: #28a745; stroke: #1e7e34; }}';
            svg += '.gene-without-guides {{ fill: #dc3545; stroke: #bd2130; }}';
            svg += '.gene-unknown {{ fill: #6c757d; stroke: #545b62; }}';
            svg += '.metabolite-node {{ fill: #17a2b8; stroke: #117a8b; stroke-width: 2; }}';
            svg += '.interaction-line {{ stroke: #333; stroke-width: 1.5; fill: none; }}';
            svg += '.node-label {{ font-family: Arial, sans-serif; font-size: 12px; text-anchor: middle; pointer-events: none; }}';
            svg += '</style></defs>';
            
            // 代謝物ノード（四角形）をいくつか追加
            const metabolites = [
                {{ name: 'PRPP', x: 100, y: 100 }},
                {{ name: 'ATP', x: 150, y: 120 }},
                {{ name: 'Histidine', x: 500, y: 450 }}
            ];
            
            metabolites.forEach(met => {{
                svg += `<rect x="${{met.x - 25}}" y="${{met.y - 10}}" width="50" height="20" class="metabolite-node" title="${{met.name}}" />`;
                svg += `<text x="${{met.x}}" y="${{met.y + 3}}" class="node-label" fill="#000">${{met.name}}</text>`;
            }});
            
            // 相互作用線をいくつか追加
            for (let i = 0; i < genes.length - 1; i++) {{
                svg += `<line x1="${{genes[i].x}}" y1="${{genes[i].y}}" x2="${{genes[i+1].x}}" y2="${{genes[i+1].y}}" class="interaction-line" />`;
            }}
            
            // 遺伝子ノード（楕円形）を追加
            genes.forEach(gene => {{
                const nodeClass = gene.hasGuides ? 'gene-with-guides' : 'gene-without-guides';
                svg += `<ellipse cx="${{gene.x}}" cy="${{gene.y}}" rx="30" ry="15" class="gene-node ${{nodeClass}}" `;
                svg += `data-gene-database="SGD" data-gene-id="${{gene.sgd_id}}" `;
                svg += `data-gene-key="SGD::${{gene.sgd_id}}" data-sgd-id="${{gene.sgd_id}}" `;
                svg += `data-gene-name="${{gene.name}}" title="${{gene.name}} (SGD:${{gene.sgd_id}})" />`;
                svg += `<text x="${{gene.x}}" y="${{gene.y + 3}}" class="node-label" fill="#000">${{gene.name}}</text>`;
            }});
            
            // 凡例を追加
            const legendY = height - 80;
            svg += '<g id="legend">';
            svg += `<rect x="10" y="${{legendY}}" width="200" height="70" fill="white" stroke="#ccc" stroke-width="1" opacity="0.9" />`;
            svg += `<text x="20" y="${{legendY + 15}}" class="node-label" font-weight="bold">Legend</text>`;
            svg += `<ellipse cx="30" cy="${{legendY + 30}}" rx="8" ry="6" class="gene-node gene-with-guides" />`;
            svg += `<text x="45" y="${{legendY + 33}}" class="node-label" font-size="9px">With Guides</text>`;
            svg += `<ellipse cx="30" cy="${{legendY + 45}}" rx="8" ry="6" class="gene-node gene-without-guides" />`;
            svg += `<text x="45" y="${{legendY + 48}}" class="node-label" font-size="9px">Without Guides</text>`;
            svg += `<rect x="120" y="${{legendY + 39}}" width="16" height="12" class="metabolite-node" />`;
            svg += `<text x="140" y="${{legendY + 48}}" class="node-label" font-size="9px">Metabolite</text>`;
            svg += '</g>';
            
            svg += '</svg>';
            return svg;
        }}
        
        function addGeneClickEvents() {{
            const geneNodes = document.querySelectorAll('.gene-node');
            geneNodes.forEach(node => {{
                node.addEventListener('click', function() {{
                    const geneName = this.getAttribute('data-gene-name');
                    const geneDatabase = this.getAttribute('data-gene-database');
                    const geneId = this.getAttribute('data-gene-id');
                    const geneKey = this.getAttribute('data-gene-key');
                    const sgdId = this.getAttribute('data-sgd-id');
                    const hasGuides = this.classList.contains('gene-with-guides');
                    
                    // 詳細なモーダルを表示
                    showGeneDetails(geneName, geneDatabase, geneId, geneKey, sgdId, hasGuides);
                }});
            }});
        }}
        
        function showGeneDetails(geneName, geneDatabase, geneId, geneKey, sgdId, hasGuides) {{
            try {{
                // 埋め込みデータからパスウェイガイドデータを取得
                const dataElement = document.getElementById('pathway-data');
                const data = JSON.parse(dataElement.textContent);
                
                // 現在表示中のパスウェイのガイドデータを取得
                const pathwayGuideData = data.pathway_guide_data[currentPathwayId];
                const pathwaySvgData = data.pathway_svgs[currentPathwayId];
                const geneStatus = pathwaySvgData && pathwaySvgData.guide_status_map
                    ? (pathwaySvgData.guide_status_map[geneKey] || {{}})
                    : {{}};
                
                // モーダルタイトルを設定
                document.getElementById('gene-detail-title').textContent = `${{geneName}} - Gene Details`;
                
                // 遺伝子基本情報を表示
                const basicInfoDiv = document.getElementById('gene-basic-info');
                basicInfoDiv.innerHTML = `
                    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem;">
                        <div>
                            <strong>Gene Name:</strong> ${{geneName || 'N/A'}}
                        </div>
                        <div>
                            <strong>Pathway Xref Database:</strong> ${{geneDatabase || 'Unavailable'}}
                        </div>
                        <div>
                            <strong>Pathway Xref ID:</strong> ${{geneId || 'Unavailable'}}
                        </div>
                        <div>
                            <strong>PW-KFI Gene Key:</strong> ${{geneStatus.pw_kfi_gene_key || 'Unavailable'}}
                        </div>
                        <div>
                            <strong>Source Xrefs:</strong> ${{(geneStatus.source_xref_keys || []).join(', ') || 'Unavailable'}}
                        </div>
                        <div>
                            <strong>Resolved SGD annotation ID:</strong> ${{sgdId || 'Unavailable'}}
                        </div>
                        <div>
                            <strong>Match Status:</strong> ${{geneStatus.match_status || 'unavailable'}}
                        </div>
                        <div>
                            <strong>Failure Reason:</strong> ${{geneStatus.failure_reason || '—'}}
                        </div>
                        <div>
                            <strong>Guide Status:</strong> 
                            <span style="color: ${{hasGuides ? '#28a745' : '#dc3545'}}; font-weight: bold;">
                                ${{hasGuides ? '✓ Available' : '✗ Not Available'}}
                            </span>
                        </div>
                    </div>
                `;
                
                // ガイドRNA詳細テーブルを構築
                const tbody = document.getElementById('gene-guides-table-body');
                const emptyDiv = document.getElementById('gene-guides-empty');
                const tableContainer = document.querySelector('#gene-guides-table').parentElement;
                
                tbody.innerHTML = '';
                
                // 該当遺伝子のガイドRNAを検索
                let geneGuides = [];
                if (pathwayGuideData && pathwayGuideData.guides) {{
                    geneGuides = pathwayGuideData.guides.filter(
                        guide => guide.gene_key === geneKey ||
                            (Array.isArray(guide.source_xref_keys) &&
                                guide.source_xref_keys.includes(geneKey))
                    );
                }}
                
                if (geneGuides.length === 0) {{
                    // ガイドRNAがない場合
                    tableContainer.style.display = 'none';
                    emptyDiv.style.display = 'block';
                }} else {{
                    // ガイドRNA詳細を表示
                    tableContainer.style.display = 'block';
                    emptyDiv.style.display = 'none';
                    
                    geneGuides.forEach(guide => {{
                        const row = tbody.insertRow();
                        row.innerHTML = `
                            <td style="font-family: monospace; font-size: 0.9em; background: #f8f9fa; padding: 0.3rem; cursor: pointer;"
                                title="Click to copy: ${{guide.guide_sequence || ''}}"
                                onclick="copySequenceToClipboard(this, '${{guide.guide_sequence || ''}}')"
                            >${{guide.guide_sequence || 'N/A'}}</td>
                            <td>${{guide.start_pos || 'N/A'}}</td>
                            <td>${{guide.end_pos || 'N/A'}}</td>
                            <td>${{guide['MM0'] || 'N/A'}}</td>
                            <td>${{guide['MM1-2'] || 'N/A'}}</td>
                            <td>${{guide['off_target_match_score'] || 'N/A'}}</td>
                            <td>${{guide['strand'] || 'N/A'}}</td>
                        `;
                    }});
                    
                    // 統計情報を追加
                    const statsRow = tbody.insertRow();
                    statsRow.style.backgroundColor = '#e9ecef';
                    statsRow.style.fontWeight = 'bold';
                    statsRow.innerHTML = `
                        <td colspan="6" style="text-align: center; padding: 1rem;">
                            <strong>統計:</strong> ${{geneGuides.length}} Guide RNA design completed
                        </td>
                    `;
                }}
                
                // モーダルを表示
                document.getElementById('gene-detail-modal').style.display = 'block';
                
            }} catch (error) {{
                console.error('Genetic detail display error:', error);
                alert('Failed to display detailed gene information.: ' + error.message);
            }}
        }}
        
        function copySequenceToClipboard(element, sequence) {{
            if (sequence && sequence !== 'N/A') {{
                navigator.clipboard.writeText(sequence).then(() => {{
                    // 視覚的フィードバック
                    const originalBg = element.style.backgroundColor;
                    const originalColor = element.style.color;
                    element.style.backgroundColor = '#28a745';
                    element.style.color = 'white';
                    
                    setTimeout(() => {{
                        element.style.backgroundColor = originalBg;
                        element.style.color = originalColor;
                    }}, 1000);
                    
                    console.log('Sequence copied:', sequence);
                }}).catch(err => {{
                    console.error('Clipboard copy failed:', err);
                    // フォールバック処理
                    const textArea = document.createElement('textarea');
                    textArea.value = sequence;
                    document.body.appendChild(textArea);
                    textArea.select();
                    document.execCommand('copy');
                    document.body.removeChild(textArea);
                    
                    // 視覚的フィードバック
                    const originalBg = element.style.backgroundColor;
                    element.style.backgroundColor = '#28a745';
                    element.style.color = 'white';
                    setTimeout(() => {{
                        element.style.backgroundColor = originalBg;
                        element.style.color = '';
                    }}, 1000);
                }});
            }}
        }}
        
        function zoomDiagram(factor) {{
            if (!currentPathwayId) return;
            
            currentZoom *= factor;
            // ズーム範囲を制限
            currentZoom = Math.max(0.2, Math.min(5.0, currentZoom));
            
            updateDiagramZoom();
        }}
        
        function resetZoom() {{
            if (!currentPathwayId) return;
            
            currentZoom = 1.0;
            updateDiagramZoom();
        }}
        
        function updateDiagramZoom() {{
            const svgElement = document.querySelector('#pathway-diagram-content svg');
            if (svgElement) {{
                svgElement.style.transform = `scale(${{currentZoom}})`;
                svgElement.style.transformOrigin = '0 0';
                
                // コンテナのサイズを調整
                const container = document.getElementById('pathway-diagram-content');
                const originalWidth = parseFloat(svgElement.getAttribute('width')) || 600;
                const originalHeight = parseFloat(svgElement.getAttribute('height')) || 500;
                
                container.style.width = (originalWidth * currentZoom) + 'px';
                container.style.height = (originalHeight * currentZoom) + 'px';
            }}
        }}
        
        function downloadPathwayCSV() {{
            // CSV ダウンロード機能
            const csvContent = generateCSVContent();
            const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
            const link = document.createElement('a');
            const url = URL.createObjectURL(blob);
            link.setAttribute('href', url);
            link.setAttribute('download', 'pathway_analysis.csv');
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }}
        
        function generateCSVContent() {{
            const headers = [
                'Pathway ID', 'Pathway Name', 'Raw GPML Nodes',
                'Evaluable Nodes', 'Unique Xrefs', 'Unique PW-KFI Genes',
                'PW-KFI Genes with Guides', 'PW-KFI Score', 'Targetability',
                'Annotation Match Coverage'
            ];
            const rows = filteredData.map(pathway => [
                pathway.pathway_id,
                `"${{pathway.pathway_name}}"`,
                pathway.total_gene_nodes,
                pathway.evaluable_gene_nodes,
                pathway.unique_xref_count,
                pathway.total_genes,
                pathway.genes_with_guides,
                pathway.pw_kfi_score.toFixed(4),
                pathway.targetable_type,
                (pathway.coverage_rate * 100).toFixed(2) + '%'
            ]);
            
            return [headers.join(','), ...rows.map(row => row.join(','))].join('\\n');
        }}
        
        function downloadPathwayGuides(pathwayId) {{
            try {{
                // 埋め込みデータからパスウェイガイドデータを取得
                const dataElement = document.getElementById('pathway-data');
                const data = JSON.parse(dataElement.textContent);
                const pathwayGuideData = data.pathway_guide_data[pathwayId];
                
                if (!pathwayGuideData || !pathwayGuideData.guides || pathwayGuideData.guides.length === 0) {{
                    alert('Guide RNA information for this pathway could not be found.: ' + pathwayId);
                    return;
                }}
                
                // パスウェイガイドCSVコンテンツを生成
                const csvContent = generatePathwayGuidesCSVContent(pathwayGuideData);
                
                // CSVファイルダウンロード
                const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
                const link = document.createElement('a');
                const url = URL.createObjectURL(blob);
                
                // ファイル名を生成（パスウェイIDを含む）
                const fileName = `pathway_${{pathwayId}}_guides.csv`;
                
                link.setAttribute('href', url);
                link.setAttribute('download', fileName);
                link.style.visibility = 'hidden';
                document.body.appendChild(link);
                link.click();
                document.body.removeChild(link);
                
                // 成功通知
                console.log(`Guide RNA information for pathway ${{pathwayId}} has been downloaded.: ${{fileName}}`);
                
            }} catch (error) {{
                console.error('Guide RNA download error:', error);
                alert('Guide RNA download failed.: ' + error.message);
            }}
        }}
        
        function generatePathwayGuidesCSVContent(pathwayGuideData) {{
            // CSVヘッダー
            const headers = [
                'Pathway ID',
                'Pathway Name',
                'Pathway Xref Database',
                'Pathway Xref ID',
                'Gene Key',
                'PW-KFI Gene Key',
                'Source Xref Keys',
                'Source Databases',
                'Annotation Source Database',
                'Matched Annotation Key',
                'Matched Annotation ID',
                'Resolved Annotation Subject ID',
                'Match Method',
                'SGD ID',
                'Ensembl ID',
                'Gene Name',
                'Unique ID',
                'Description',
                'Guide Unique ID',
                'Start Position',
                'End Position',
                'Guide Sequence',
                'MM0',
                'MM1-2',
                'off_target_match_score',
                'strand',
            ];
            
            // データ行を生成
            const rows = pathwayGuideData.guides.map(guide => [
                guide.pathway_id || '',
                `"${{(guide.pathway_name || '').replace(/"/g, '""')}}"`, // CSVエスケープ
                guide.pathway_gene_database || '',
                guide.pathway_gene_id || '',
                guide.gene_key || '',
                guide.pw_kfi_gene_key || '',
                `"${{(guide.source_xref_keys || []).join(';').replace(/"/g, '""')}}"`,
                `"${{(guide.source_databases || []).join(';').replace(/"/g, '""')}}"`,
                guide.annotation_source_database || '',
                guide.matched_annotation_key || '',
                guide.matched_annotation_id || '',
                guide.resolved_annotation_subject_id || '',
                guide.match_method || '',
                guide.sgd_id || '',
                guide.ensembl_id || '',
                `"${{(guide.gene_name || '').replace(/"/g, '""')}}"`, // CSVエスケープ
                guide.unique_id || '',
                `"${{(guide.description || '').replace(/"/g, '""')}}"`, // CSVエスケープ
                guide.guide_unique_id || '',
                guide.start_pos || '',
                guide.end_pos || '',
                guide.guide_sequence || '',
                guide['MM0'] == null ? '' : guide['MM0'],
                guide['MM1-2'] == null ? '' : guide['MM1-2'],
                guide.off_target_match_score || '',
                guide.strand || '',
            ]);
            
            // CSV形式に変換
            return [headers.join(','), ...rows.map(row => row.join(','))].join('\\n');
        }}
    </script>
</body>
</html>'''
        
        return html_template
    
    def generate_pathway_html_report(self, pathways, pathway_genes, summary_data):
        """Generate a pathway HTML report"""
        html_content = self.create_pathway_html_template(pathways, summary_data, pathway_genes)
        
        # ロゴファイルのコピー
        img_dir = os.path.join(self.output_dir, 'static', 'img')
        Path(img_dir).mkdir(parents=True, exist_ok=True)
        src_logo = os.path.join('src', 'web', 'static', 'img', 'logo.png')
        dst_logo = os.path.join(img_dir, 'logo.png')
        if os.path.exists(src_logo):
            shutil.copy2(src_logo, dst_logo)
            self.logger.info(f"Logo copied: {src_logo} -> {dst_logo}")
        
        # HTMLファイル保存
        report_name = self.pathway_report_config.get('output_name', 'pathway_report.html')
        report_path = os.path.join(self.output_dir, report_name)
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        self.logger.info(f"Pathway HTML Report Generation: {report_path}")
        return report_path
    
    def run(self):
        """Full execution of pathway report generation"""
        try:
            self.logger.info("Started generating pathway report....")
            
            # データロード
            pathways, pathway_genes = self.load_pathway_analysis_results()
            pathway_gene_audit = self.load_pathway_gene_match_audit()
            
            # パスウェイ別ガイドCSVデータ生成
            pathway_guide_data = self.generate_pathway_guide_csv_data(pathways)
            
            # データベーステーブル作成・データ挿入
            db_path = self.create_pathway_database_tables()
            self.populate_pathway_database(
                pathways, pathway_genes, db_path, pathway_gene_audit
            )
            
            # JSONサマリー生成（ガイドデータとSVGデータを含む）
            summary_data = self.generate_pathway_json_summary(
                pathways, pathway_guide_data
            )
            
            # HTMLレポート生成
            report_path = self.generate_pathway_html_report(pathways, pathway_genes, summary_data)
            
            # 結果まとめ
            results = {
                'database_path': db_path,
                'summary_json_path': os.path.join(self.web_data_dir, 'pathway_summary.json'),
                'html_report_path': report_path,
                'statistics': summary_data['statistics']['pathway_analysis'],
                'pathway_guide_data': pathway_guide_data
            }
            
            self.logger.info(f"Pathway report generation complete.: {report_path}")
            return results
            
        except Exception as e:
            self.logger.error(f"An error occurred during pathway report generation.: {e}")
            raise
