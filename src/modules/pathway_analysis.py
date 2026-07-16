#!/usr/bin/env python3
"""
Pathway Analysis Module
Pathway Knock-out Feasibility Index (PW-KFI) calculation using WikiPathways GPML files
"""

import os
import logging
import xml.etree.ElementTree as ET
import re
import csv
import gzip
import sys
from collections import Counter
from pathlib import Path


try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(2 ** 31 - 1)


class PathwayAnalysis:
    """Pathway analysis processing class"""
    
    def __init__(self, config):
        """
        Initialize
        
        Args:
            config (dict): Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        self.pathway_config = self.config.get('pathway_analysis', {})
        
        # パスウェイ統計情報
        self.coverage_stats = {
            'total_pathways': 0,
            'processable_pathways': 0,
            'skipped_pathways': 0,
            'coverage_by_database': {}
        }
        
        # パスウェイデータストレージ
        self.pathway_data = {}
        self.pathway_gene_audit = []
        self.annotation_identifier_issues = []

    def _target_database_configs(self):
        """Return target database configurations keyed by exact GPML database name."""
        return {
            item['database_name']: item
            for item in self.pathway_config.get('target_gene_databases', [])
        }
        
    def validate_configuration(self):
        """Validate configuration"""
        required_keys = ['wikipathways_gpml_dir', 'target_gene_databases']
        
        for key in required_keys:
            if key not in self.pathway_config:
                raise ValueError(f"The key '{key}' required for the pathway_analysis configuration could not be found.")
        
        gpml_dir = self.pathway_config['wikipathways_gpml_dir']
        if not os.path.exists(gpml_dir):
            raise FileNotFoundError(f"GPML directory not found.: {gpml_dir}")

        database_names = set()
        for db_config in self.pathway_config['target_gene_databases']:
            for key in ('database_name', 'id_pattern', 'annotation_key'):
                if not db_config.get(key):
                    raise ValueError(
                        f"target_gene_databases entry is missing '{key}': {db_config}"
                    )

            database_name = db_config['database_name']
            if database_name in database_names:
                raise ValueError(
                    f"Duplicate target database configuration: {database_name}"
                )
            database_names.add(database_name)

            try:
                re.compile(db_config['id_pattern'])
            except re.error as exc:
                raise ValueError(
                    f"Invalid ID pattern for {database_name}: {exc}"
                ) from exc
        
        self.logger.info("Pathway analysis Validate configuration complete.")
        return True
        
    def get_gpml_files(self):
        """Get a list of GPML files"""
        gpml_dir = self.pathway_config['wikipathways_gpml_dir']
        gpml_files = list(Path(gpml_dir).glob("*.gpml"))
        
        self.logger.info(f"Number of GPML files discovered: {len(gpml_files)}")
        return gpml_files
    
    def parse_gpml_file(self, gpml_file):
        """Parse individual GPML file"""
        try:
            tree = ET.parse(gpml_file)
            root = tree.getroot()
            
            # パスウェイメタデータ抽出
            pathway_info = self.extract_pathway_metadata(root, gpml_file)
            
            # 遺伝子情報抽出
            pathway_genes, rejected_genes = self.extract_pathway_genes(
                root, include_rejected=True
            )
            
            # 統計更新
            self.coverage_stats['total_pathways'] += 1
            
            if pathway_genes:
                pathway_info['genes'] = pathway_genes
                pathway_info['evaluable_gene_node_count'] = len(pathway_genes)
                pathway_info['unique_xref_count'] = len({
                    (gene['database'], gene['gene_id'])
                    for gene in pathway_genes
                })
                # Cross-database aliases are resolved after annotation loading.
                pathway_info['total_genes'] = pathway_info['unique_xref_count']
                pathway_info['total_gene_nodes'] = len(pathway_genes) + len(rejected_genes)
                pathway_info['unmatchable_genes'] = rejected_genes
                pathway_info['unmatchable_gene_count'] = len(rejected_genes)
                self.coverage_stats['processable_pathways'] += 1
                return pathway_info
            else:
                pathway_info['genes'] = []
                pathway_info['evaluable_gene_node_count'] = 0
                pathway_info['unique_xref_count'] = 0
                pathway_info['total_genes'] = 0
                pathway_info['total_gene_nodes'] = len(rejected_genes)
                pathway_info['unmatchable_genes'] = rejected_genes
                pathway_info['unmatchable_gene_count'] = len(rejected_genes)
                self.logger.warning(
                    f"No evaluable configured-database gene was found.: {gpml_file.name}"
                )
                self.coverage_stats['skipped_pathways'] += 1
                return pathway_info
                
        except ET.ParseError as e:
            self.logger.error(f"GPML file parsing error {gpml_file}: {e}")
            self.coverage_stats['skipped_pathways'] += 1
            return None
    
    def extract_pathway_metadata(self, root, gpml_file):
        """Extract pathway metadata"""
        # XMLネームスペース処理
        namespace = {'gpml': 'http://pathvisio.org/GPML/2013a'}
        
        # パスウェイ基本情報
        pathway_name = root.get('Name', 'Unknown')
        pathway_organism = root.get('Organism', 'Unknown')
        
        # WikiPathways ID抽出
        version_info = root.get('Version', '')
        pathway_id_match = re.search(r'(WP\d+)', version_info)
        pathway_id = pathway_id_match.group(1) if pathway_id_match else os.path.basename(gpml_file).split('_')[1]
        
        # 説明文抽出
        description = ""
        comment_elements = root.findall('.//gpml:Comment[@Source="WikiPathways-description"]', namespace)
        if comment_elements:
            description = comment_elements[0].text or ""
        
        return {
            'pathway_id': pathway_id,
            'pathway_name': pathway_name,
            'pathway_organism': pathway_organism,
            'pathway_description': description,
            'gpml_file_path': str(gpml_file)
        }
    
    def extract_pathway_genes(self, root, include_rejected=False):
        """Extract strictly valid pathway genes and optionally rejected nodes."""
        namespace = {'gpml': 'http://pathvisio.org/GPML/2013a'}
        genes = []
        rejected = []
        database_configs = self._target_database_configs()
        
        # GeneProduct タイプのDataNodeを検索
        gene_nodes = root.findall('.//gpml:DataNode[@Type="GeneProduct"]', namespace)
        
        for gene_node in gene_nodes:
            xref_element = gene_node.find('gpml:Xref', namespace)
            gene_name = gene_node.get('TextLabel', 'Unknown')

            if xref_element is None:
                rejected.append({
                    'graph_id': gene_node.get('GraphId', ''),
                    'database': '',
                    'gene_id': '',
                    'gene_name': gene_name,
                    'annotation_key': '',
                    'reason': 'missing_xref',
                })
                continue

            database = xref_element.get('Database', '').strip()
            gene_id = xref_element.get('ID', '').strip()
            db_config = database_configs.get(database)

            if db_config is None:
                rejected.append({
                    'graph_id': gene_node.get('GraphId', ''),
                    'database': database,
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'annotation_key': '',
                    'reason': 'missing_database' if not database else 'unsupported_database',
                })
                continue

            if not gene_id:
                rejected.append({
                    'graph_id': gene_node.get('GraphId', ''),
                    'database': database,
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'annotation_key': db_config['annotation_key'],
                    'reason': 'missing_identifier',
                })
                continue

            if re.fullmatch(db_config['id_pattern'], gene_id) is None:
                rejected.append({
                    'graph_id': gene_node.get('GraphId', ''),
                    'database': database,
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'annotation_key': db_config['annotation_key'],
                    'reason': 'invalid_identifier_format',
                })
                continue

            genes.append({
                'graph_id': gene_node.get('GraphId', ''),
                'database': database,
                'gene_id': gene_id,
                'gene_name': gene_name,
                'annotation_key': db_config['annotation_key'],
            })

        if include_rejected:
            return genes, rejected
        return genes
    
    @staticmethod
    def _guide_entries(guide_targets):
        """Return non-empty guide entries in deterministic de-duplicated order."""
        return sorted(set(
            item.strip()
            for item in (guide_targets or '').split(',')
            if item.strip()
        ))

    def _load_sgd_identifier_crosswalk(self):
        """Load an exact Ensembl-systematic-ID to SGD subject-ID cross-reference."""
        resolver = {
            'available': False,
            'source_path': '',
            'ensembl_to_subject': {},
            'ambiguous_ensembl_ids': set(),
            'failure_reason': 'resolver_not_configured',
        }

        ensembl_config = self._target_database_configs().get('Ensembl')
        if ensembl_config is None:
            return resolver

        sgd_path = self.config.get('input', {}).get('sgd_db', '')
        if not sgd_path:
            resolver['failure_reason'] = 'resolver_source_not_configured'
            return resolver

        resolver['source_path'] = os.path.abspath(sgd_path)
        if not os.path.isfile(sgd_path):
            resolver['failure_reason'] = 'resolver_source_unavailable'
            return resolver

        opener = gzip.open if str(sgd_path).endswith('.gz') else open
        mapping = {}
        ambiguous = set()
        header_pattern = re.compile(
            r'^>(?P<ensembl_id>[^|\s]+)\|[^|\s]*\|'
            r'(?P<subject_id>SGDID:S\d{9})(?:\s|$)'
        )

        try:
            with opener(sgd_path, 'rt', encoding='utf-8') as handle:
                for line_number, line in enumerate(handle, 1):
                    if not line.startswith('>'):
                        continue

                    match = header_pattern.match(line.rstrip('\r\n'))
                    if match is None:
                        continue

                    ensembl_id = match.group('ensembl_id')
                    subject_id = match.group('subject_id')
                    if re.fullmatch(ensembl_config['id_pattern'], ensembl_id) is None:
                        continue

                    previous = mapping.get(ensembl_id)
                    if previous is None:
                        mapping[ensembl_id] = subject_id
                    elif previous != subject_id:
                        ambiguous.add(ensembl_id)
                        mapping.pop(ensembl_id, None)
                        self.annotation_identifier_issues.append({
                            'row_number': line_number,
                            'database': 'Ensembl',
                            'annotation_key': 'ensembl_id',
                            'identifier': ensembl_id,
                            'reason': 'ambiguous_cross_reference',
                            'query_id': '',
                            'unique_id': '',
                            'source_path': resolver['source_path'],
                        })
        except (OSError, UnicodeError) as exc:
            resolver['failure_reason'] = 'resolver_source_unreadable'
            self.logger.warning(
                f"Unable to load exact SGD identifier cross-reference: {exc}"
            )
            return resolver

        if not mapping and not ambiguous:
            resolver['failure_reason'] = 'resolver_contains_no_valid_cross_references'
            self.logger.warning(
                "The configured SGD database contains no valid exact Ensembl-to-SGD "
                "cross-references for the configured identifier pattern."
            )
            return resolver

        resolver.update({
            'available': True,
            'ensembl_to_subject': mapping,
            'ambiguous_ensembl_ids': ambiguous,
            'failure_reason': '',
        })
        self.logger.info(
            "Exact SGD identifier cross-reference loaded: "
            f"{len(mapping)} Ensembl systematic IDs, {len(ambiguous)} ambiguous"
        )
        return resolver

    def load_annotation_results(self):
        """Load annotation results into strict database-specific indexes."""
        annotation_file = os.path.join(self.intermediate_dir, 'annotation_with_guides_final.tsv')
        
        if not os.path.exists(annotation_file):
            raise FileNotFoundError(f"Annotation result file not found.: {annotation_file}")

        self.annotation_identifier_issues = []
        database_configs = self._target_database_configs()
        annotation_data = {
            'by_database': {name: {} for name in database_configs},
            'ambiguous_by_database': {name: set() for name in database_configs},
            'database_status': {},
            'by_subject_id': {},
            'ambiguous_subject_ids': set(),
            'subject_id_column_available': False,
            'cross_references': {},
        }

        candidate_records = []
        issue_keys = set()

        def add_issue(issue):
            """Append one deterministic issue per semantic conflict."""
            issue = {
                'row_number': '',
                'database': '',
                'annotation_key': '',
                'identifier': '',
                'reason': '',
                'query_id': '',
                'unique_id': '',
                'conflicting_fields': '',
                'candidate_values': '',
                'source_path': os.path.abspath(annotation_file),
                **issue,
            }
            key = tuple(
                str(issue.get(field, ''))
                for field in (
                    'source_path', 'database', 'annotation_key', 'identifier',
                    'reason', 'query_id', 'unique_id', 'conflicting_fields',
                    'candidate_values',
                )
            )
            if key not in issue_keys:
                issue_keys.add(key)
                self.annotation_identifier_issues.append(issue)

        with open(annotation_file, 'r', encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            fieldnames = set(reader.fieldnames or [])
            annotation_data['subject_id_column_available'] = 'subject_id' in fieldnames

            for database, db_config in database_configs.items():
                annotation_key = db_config['annotation_key']
                column_available = annotation_key in fieldnames
                annotation_data['database_status'][database] = {
                    'annotation_key': annotation_key,
                    'column_available': column_available,
                }
                if not column_available:
                    add_issue({
                        'database': database,
                        'annotation_key': annotation_key,
                        'reason': 'annotation_identifier_column_unavailable',
                    })

            for row_number, row in enumerate(reader, 2):
                unique_id = (row.get('unique_id') or '').strip()
                subject_id = (row.get('subject_id') or '').strip()
                if not unique_id or not subject_id:
                    add_issue({
                        'row_number': row_number,
                        'annotation_key': 'subject_id',
                        'identifier': subject_id,
                        'reason': 'missing_unique_id' if not unique_id else 'missing_subject_id',
                        'query_id': (row.get('query_id') or '').strip(),
                        'unique_id': unique_id,
                    })
                    continue

                guides = self._guide_entries(row.get('guide_targets', ''))
                annotation_source_database = (
                    row.get('annotation_source_database') or ''
                ).strip()
                if (
                    not annotation_source_database
                    and re.fullmatch(r'SGDID:S\d{9}', subject_id)
                ):
                    annotation_source_database = 'SGD'
                if not annotation_source_database:
                    annotation_source_database = 'Unspecified'
                record = dict(row)
                record.update({
                    'unique_id': unique_id,
                    'subject_id': subject_id,
                    'gene_name': (row.get('gene_name') or '').strip(),
                    'description': (row.get('description') or '').strip(),
                    'annotation_source_database': annotation_source_database,
                    'guide_targets': ','.join(guides),
                    'guide_count': len(guides),
                    'has_guides': bool(guides),
                    'annotation_query_ids': [
                        (row.get('query_id') or '').strip()
                    ] if (row.get('query_id') or '').strip() else [],
                    '_annotation_row_number': row_number,
                })

                for database, db_config in database_configs.items():
                    annotation_key = db_config['annotation_key']
                    if annotation_key not in fieldnames:
                        continue

                    identifier = (row.get(annotation_key) or '').strip()
                    if not identifier:
                        continue
                    if re.fullmatch(db_config['id_pattern'], identifier) is None:
                        add_issue({
                            'row_number': row_number,
                            'database': database,
                            'annotation_key': annotation_key,
                            'identifier': identifier,
                            'reason': 'invalid_annotation_identifier_format',
                            'query_id': (row.get('query_id') or '').strip(),
                            'unique_id': unique_id,
                        })
                candidate_records.append(record)

        records_by_subject = {}
        for record in candidate_records:
            records_by_subject.setdefault(record['subject_id'], []).append(record)

        identity_fields = sorted({
            'unique_id', 'annotation_source_database', 'classification',
            'sgd_id', 'ensembl_id', 'gene_name', 'description',
            *(config['annotation_key'] for config in database_configs.values()),
        })

        def nonempty_values(records, field):
            return sorted({
                str(record.get(field, '') or '').strip()
                for record in records
                if str(record.get(field, '') or '').strip()
            })

        def numeric_value(value, default):
            try:
                number = float(str(value).strip())
            except (TypeError, ValueError):
                return default
            if number != number or number in (float('inf'), float('-inf')):
                return default
            return number

        def representative_key(record):
            """Rank query-specific alignments independently of input row order."""
            return (
                numeric_value(record.get('evalue'), float('inf')),
                -numeric_value(record.get('bit_score'), float('-inf')),
                -numeric_value(record.get('identity'), float('-inf')),
                -numeric_value(record.get('length'), float('-inf')),
                str(record.get('query_id', '') or ''),
                tuple(
                    (key, str(record.get(key, '') or ''))
                    for key in sorted(record)
                    if key != '_annotation_row_number'
                ),
            )

        for subject_id in sorted(records_by_subject):
            subject_records = records_by_subject[subject_id]
            conflicts = {
                field: nonempty_values(subject_records, field)
                for field in identity_fields
                if len(nonempty_values(subject_records, field)) > 1
            }
            query_ids = sorted({
                str(record.get('query_id', '') or '').strip()
                for record in subject_records
                if str(record.get('query_id', '') or '').strip()
            })
            unique_ids = nonempty_values(subject_records, 'unique_id')
            if conflicts:
                annotation_data['ambiguous_subject_ids'].add(subject_id)
                for database, db_config in database_configs.items():
                    annotation_key = db_config['annotation_key']
                    for identifier in nonempty_values(
                        subject_records, annotation_key
                    ):
                        if re.fullmatch(
                            db_config['id_pattern'], identifier
                        ) is not None:
                            annotation_data['ambiguous_by_database'][
                                database
                            ].add(identifier)
                add_issue({
                    'database': 'Annotation',
                    'annotation_key': 'subject_id',
                    'identifier': subject_id,
                    'reason': 'ambiguous_annotation_subject_id',
                    'query_id': ','.join(query_ids),
                    'unique_id': ','.join(unique_ids),
                    'conflicting_fields': ','.join(sorted(conflicts)),
                    'candidate_values': ';'.join(
                        f"{field}={'|'.join(conflicts[field])}"
                        for field in sorted(conflicts)
                    ),
                })
                continue

            aggregate = min(subject_records, key=representative_key).copy()
            for field in identity_fields:
                values = nonempty_values(subject_records, field)
                aggregate[field] = values[0] if values else ''
            guides = sorted({
                guide
                for record in subject_records
                for guide in self._guide_entries(record.get('guide_targets', ''))
            })
            aggregate.update({
                'subject_id': subject_id,
                'query_id': str(aggregate.get('query_id', '') or '').strip(),
                'annotation_query_ids': query_ids,
                'annotation_row_count': len(subject_records),
                'guide_targets': ','.join(guides),
                'guide_count': len(guides),
                'has_guides': bool(guides),
            })
            aggregate.pop('_annotation_row_number', None)
            explicit_sgd_id = str(aggregate.get('sgd_id', '') or '')
            subject_sgd_match = re.fullmatch(
                r'SGDID:(S\d{9})', subject_id
            )
            if (
                subject_sgd_match is not None
                and explicit_sgd_id
                and explicit_sgd_id != subject_sgd_match.group(1)
            ):
                annotation_data['ambiguous_subject_ids'].add(subject_id)
                for database, db_config in database_configs.items():
                    annotation_key = db_config['annotation_key']
                    identifier = str(
                        aggregate.get(annotation_key, '') or ''
                    )
                    if identifier and re.fullmatch(
                        db_config['id_pattern'], identifier
                    ) is not None:
                        annotation_data['ambiguous_by_database'][database].add(
                            identifier
                        )
                add_issue({
                    'database': 'Annotation',
                    'annotation_key': 'subject_id,sgd_id',
                    'identifier': subject_id,
                    'reason': 'conflicting_annotation_identity',
                    'query_id': ','.join(query_ids),
                    'unique_id': ','.join(unique_ids),
                    'conflicting_fields': 'subject_id,sgd_id',
                    'candidate_values': (
                        f"subject_id={subject_id};sgd_id={explicit_sgd_id}"
                    ),
                })
                continue
            annotation_data['by_subject_id'][subject_id] = aggregate

        for database, db_config in sorted(database_configs.items()):
            annotation_key = db_config['annotation_key']
            if annotation_key not in fieldnames:
                continue

            records_by_identifier = {}
            for subject_id in sorted(annotation_data['by_subject_id']):
                record = annotation_data['by_subject_id'][subject_id]
                identifier = str(record.get(annotation_key, '') or '').strip()
                if not identifier:
                    continue
                if re.fullmatch(db_config['id_pattern'], identifier) is None:
                    continue
                records_by_identifier.setdefault(identifier, []).append(record)

            for identifier in sorted(records_by_identifier):
                records = records_by_identifier[identifier]
                identities = {
                    (
                        record.get('subject_id', ''),
                        record.get('unique_id', ''),
                        record.get('annotation_source_database', ''),
                    )
                    for record in records
                }
                if len(identities) != 1:
                    annotation_data['ambiguous_by_database'][database].add(
                        identifier
                    )
                    add_issue({
                        'database': database,
                        'annotation_key': annotation_key,
                        'identifier': identifier,
                        'reason': 'ambiguous_annotation_identifier',
                        'query_id': ','.join(sorted({
                            query_id
                            for record in records
                            for query_id in record.get(
                                'annotation_query_ids', []
                            )
                        })),
                        'unique_id': ','.join(sorted({
                            record.get('unique_id', '') for record in records
                        })),
                        'candidate_values': ';'.join(
                            '|'.join(identity) for identity in sorted(identities)
                        ),
                    })
                    continue
                annotation_data['by_database'][database][identifier] = records[0]

        annotation_data['cross_references']['Ensembl'] = (
            self._load_sgd_identifier_crosswalk()
        )
        self.annotation_identifier_issues.sort(key=lambda issue: tuple(
            str(issue.get(field, ''))
            for field in (
                'source_path', 'database', 'annotation_key', 'identifier',
                'reason', 'query_id', 'unique_id', 'conflicting_fields',
                'candidate_values', 'row_number',
            )
        ))
        self.logger.info(
            "Annotation results loading complete: "
            f"{len(annotation_data['by_subject_id'])} exact subject IDs"
        )
        return annotation_data

    @staticmethod
    def _source_gene_key(database, gene_id):
        """Return an exact GPML database-and-identifier key."""
        return f"{database}::{gene_id}"

    @staticmethod
    def _subject_gene_key(subject_id):
        """Return the canonical PW-KFI key for one resolved subject."""
        namespace = 'SGD' if subject_id.startswith('SGDID:') else 'Subject'
        return f"{namespace}::{subject_id}"

    def _canonical_subject_for_pathway_gene(self, gene, annotation_data):
        """Resolve a GPML gene to a subject only through exact evidence."""
        database = gene['database']
        gene_id = gene['gene_id']
        if database == 'SGD':
            return f"SGDID:{gene_id}", 'exact_sgd_identifier', True
        if database == 'Ensembl':
            resolver = annotation_data['cross_references'].get('Ensembl', {})
            if not resolver.get('available'):
                return '', resolver.get('failure_reason') or 'resolver_unavailable', False
            if gene_id in resolver.get('ambiguous_ensembl_ids', set()):
                return '', 'ambiguous_cross_reference', False
            subject_id = resolver.get('ensembl_to_subject', {}).get(gene_id, '')
            if not subject_id:
                return '', 'cross_reference_not_found', False
            return subject_id, 'exact_sgd_fasta_cross_reference', True
        return '', 'unsupported_canonical_namespace', False
    
    def match_pathway_genes_with_annotation(self, pathway_genes, annotation_data):
        """Match pathway genes using only exact, provenance-preserving identifiers."""
        matched_genes = []
        unmatched_genes = []
        
        for gene in pathway_genes:
            database = gene['database']
            gene_id = gene['gene_id']
            annotation_key = gene['annotation_key']
            canonical_subject_id, canonicalization_method, canonicalization_complete = (
                self._canonical_subject_for_pathway_gene(gene, annotation_data)
            )
            record = None
            match_method = ''
            mapping_source = ''
            matched_annotation_key = ''
            matched_annotation_id = ''
            failure_reason = ''

            if gene_id in annotation_data['ambiguous_by_database'].get(database, set()):
                failure_reason = 'ambiguous_annotation_identifier'
            else:
                record = annotation_data['by_database'].get(database, {}).get(gene_id)
                if record is not None:
                    match_method = 'exact_annotation_identifier'
                    mapping_source = 'annotation_with_guides_final.tsv'
                    matched_annotation_key = annotation_key
                    matched_annotation_id = gene_id

                    record_subject_id = (record.get('subject_id') or '').strip()
                    annotation_source = record.get(
                        'annotation_source_database', ''
                    )
                    if (
                        database == 'SGD'
                        and (
                            annotation_source == 'SGD'
                            or record_subject_id.startswith('SGDID:')
                        )
                        and record_subject_id != f"SGDID:{gene_id}"
                    ):
                        record = None
                        failure_reason = 'conflicting_annotation_identifier'
                    elif database == 'Ensembl':
                        resolver = annotation_data['cross_references'].get(
                            'Ensembl', {}
                        )
                        expected_subject_id = resolver.get(
                            'ensembl_to_subject', {}
                        ).get(gene_id)
                        if resolver.get('available'):
                            if gene_id in resolver.get(
                                'ambiguous_ensembl_ids', set()
                            ):
                                record = None
                                failure_reason = 'ambiguous_cross_reference'
                            elif not expected_subject_id:
                                record = None
                                failure_reason = 'cross_reference_not_found'
                            elif (
                                record.get('sgd_id')
                                and f"SGDID:{record.get('sgd_id')}"
                                != expected_subject_id
                            ):
                                record = None
                                failure_reason = 'conflicting_cross_reference'
                            elif (
                                annotation_source == 'SGD'
                                or record_subject_id.startswith('SGDID:')
                            ) and record_subject_id != expected_subject_id:
                                record = None
                                failure_reason = 'conflicting_cross_reference'

            if record is None and not failure_reason and database == 'SGD':
                subject_id = f"SGDID:{gene_id}"
                if subject_id in annotation_data['ambiguous_subject_ids']:
                    failure_reason = 'ambiguous_annotation_subject_id'
                else:
                    record = annotation_data['by_subject_id'].get(subject_id)
                    if record is not None:
                        match_method = 'exact_subject_id'
                        mapping_source = 'annotation_with_guides_final.tsv'
                        matched_annotation_key = 'subject_id'
                        matched_annotation_id = subject_id
                    elif not annotation_data['subject_id_column_available']:
                        failure_reason = 'annotation_subject_id_column_unavailable'
                    else:
                        failure_reason = 'annotation_not_found'

            if record is None and not failure_reason and database == 'Ensembl':
                resolver = annotation_data['cross_references'].get('Ensembl', {})
                if not resolver.get('available'):
                    failure_reason = resolver.get('failure_reason') or 'resolver_unavailable'
                elif gene_id in resolver.get('ambiguous_ensembl_ids', set()):
                    failure_reason = 'ambiguous_cross_reference'
                else:
                    subject_id = resolver.get('ensembl_to_subject', {}).get(gene_id)
                    if not subject_id:
                        failure_reason = 'cross_reference_not_found'
                    elif subject_id in annotation_data['ambiguous_subject_ids']:
                        failure_reason = 'ambiguous_annotation_subject_id'
                    else:
                        record = annotation_data['by_subject_id'].get(subject_id)
                        if record is None:
                            failure_reason = (
                                'annotation_not_found'
                                if annotation_data['subject_id_column_available']
                                else 'annotation_subject_id_column_unavailable'
                            )
                        else:
                            match_method = 'exact_sgd_fasta_cross_reference'
                            mapping_source = resolver.get('source_path', '')
                            matched_annotation_key = 'subject_id'
                            matched_annotation_id = subject_id

            if record is None:
                if not failure_reason:
                    database_status = annotation_data['database_status'].get(database, {})
                    if not database_status.get('column_available', False):
                        failure_reason = 'annotation_identifier_column_unavailable'
                    else:
                        failure_reason = 'no_exact_annotation_match'
                unmatched = gene.copy()
                unmatched.update({
                    'reason': failure_reason,
                    'pathway_gene_name': gene['gene_name'],
                    'pathway_gene_database': database,
                    'pathway_gene_id': gene_id,
                    'gene_key': self._source_gene_key(database, gene_id),
                    'canonical_subject_id': canonical_subject_id,
                    'pw_kfi_gene_key': (
                        self._subject_gene_key(canonical_subject_id)
                        if canonical_subject_id
                        else self._source_gene_key(database, gene_id)
                    ),
                    'canonicalization_method': canonicalization_method,
                    'canonicalization_complete': canonicalization_complete,
                    'match_status': 'unmatched',
                })
                unmatched_genes.append(unmatched)
                continue

            gene_data = record.copy()
            sgd_id = (record.get('sgd_id') or '').strip()
            if not sgd_id and matched_annotation_id.startswith('SGDID:'):
                sgd_id = matched_annotation_id[len('SGDID:'):]
            if not sgd_id:
                subject_id = (record.get('subject_id') or '').strip()
                if subject_id.startswith('SGDID:'):
                    sgd_id = subject_id[len('SGDID:'):]
            if not sgd_id and canonical_subject_id.startswith('SGDID:'):
                sgd_id = canonical_subject_id[len('SGDID:'):]
            ensembl_id = (record.get('ensembl_id') or '').strip()
            if not ensembl_id and database == 'Ensembl':
                ensembl_id = gene_id

            resolved_annotation_subject_id = record.get('subject_id', '')
            explicit_annotation_sgd_id = str(
                record.get('sgd_id', '') or ''
            ).strip()
            explicit_annotation_subject_id = (
                f"SGDID:{explicit_annotation_sgd_id}"
                if re.fullmatch(r'S\d{9}', explicit_annotation_sgd_id)
                else ''
            )
            subject_annotation_cross_reference = (
                database == 'Ensembl'
                and re.fullmatch(
                    r'SGDID:S\d{9}', str(resolved_annotation_subject_id)
                ) is not None
            )
            explicit_annotation_cross_reference = (
                database == 'Ensembl'
                and bool(explicit_annotation_subject_id)
            )
            canonical_gene_complete = (
                canonicalization_complete
                or explicit_annotation_cross_reference
                or subject_annotation_cross_reference
            )
            canonical_gene_subject_id = (
                canonical_subject_id
                if canonicalization_complete and canonical_subject_id
                else (
                    explicit_annotation_subject_id
                    or resolved_annotation_subject_id
                )
            )
            gene_data.update({
                'graph_id': gene.get('graph_id', ''),
                'pathway_gene_name': gene['gene_name'],
                'pathway_gene_database': database,
                'pathway_gene_id': gene_id,
                'gene_key': self._source_gene_key(database, gene_id),
                'annotation_key': annotation_key,
                'annotation_source_database': record.get(
                    'annotation_source_database', 'Unspecified'
                ),
                'matched_annotation_key': matched_annotation_key,
                'matched_annotation_id': matched_annotation_id,
                'resolved_annotation_subject_id': resolved_annotation_subject_id,
                'canonical_subject_id': canonical_gene_subject_id,
                'pw_kfi_gene_key': self._subject_gene_key(
                    canonical_gene_subject_id
                ),
                'canonicalization_method': (
                    canonicalization_method
                    if canonicalization_complete and canonical_subject_id
                    else (
                        'exact_annotation_sgd_id'
                        if explicit_annotation_cross_reference
                        else 'exact_annotation_subject_id'
                    )
                ),
                'canonicalization_complete': canonical_gene_complete,
                'match_method': match_method,
                'match_status': 'matched',
                'mapping_source': mapping_source,
                'sgd_id': sgd_id,
                'ensembl_id': ensembl_id,
                'original_gene_id': gene_id,
            })
            matched_genes.append(gene_data)

        return matched_genes, unmatched_genes

    @staticmethod
    def _pathway_member_sort_key(gene):
        """Return a stable order for raw GPML gene records."""
        return (
            str(gene.get('pathway_gene_database', gene.get('database', ''))),
            str(gene.get('pathway_gene_id', gene.get('gene_id', ''))),
            str(gene.get('pathway_gene_name', gene.get('gene_name', ''))),
            str(gene.get('graph_id', '')),
            str(gene.get('query_id', '')),
        )

    def collapse_pathway_gene_units(self, matched_genes, unmatched_genes):
        """Collapse GPML nodes to exact, provenance-preserving PW-KFI genes."""
        grouped_members = {}
        for status, genes in (
            ('matched', matched_genes), ('unmatched', unmatched_genes)
        ):
            for gene in genes:
                pw_kfi_gene_key = gene.get('pw_kfi_gene_key', '')
                if not pw_kfi_gene_key:
                    raise ValueError(
                        "Pathway gene is missing a canonical PW-KFI key: "
                        f"{gene!r}"
                    )
                grouped_members.setdefault(pw_kfi_gene_key, []).append(
                    (status, gene)
                )

        matched_units = []
        unmatched_units = []
        for pw_kfi_gene_key in sorted(grouped_members):
            members = grouped_members[pw_kfi_gene_key]
            members.sort(key=lambda item: self._pathway_member_sort_key(item[1]))
            matched_members = [
                gene for status, gene in members if status == 'matched'
            ]

            xref_members = {}
            for status, gene in members:
                database = gene.get(
                    'pathway_gene_database', gene.get('database', '')
                )
                gene_id = gene.get(
                    'pathway_gene_id', gene.get('gene_id', '')
                )
                source_key = self._source_gene_key(database, gene_id)
                xref_members.setdefault(source_key, []).append((status, gene))

            pathway_xrefs = []
            for source_key in sorted(xref_members):
                source_members = xref_members[source_key]
                source_members.sort(
                    key=lambda item: self._pathway_member_sort_key(item[1])
                )
                representative_status, representative = source_members[0]
                statuses = {status for status, _ in source_members}
                status = 'matched' if 'matched' in statuses else 'unmatched'
                database = representative.get(
                    'pathway_gene_database', representative.get('database', '')
                )
                gene_id = representative.get(
                    'pathway_gene_id', representative.get('gene_id', '')
                )
                names = sorted({
                    str(gene.get(
                        'pathway_gene_name', gene.get('gene_name', '')
                    ) or '')
                    for _, gene in source_members
                    if str(gene.get(
                        'pathway_gene_name', gene.get('gene_name', '')
                    ) or '')
                })
                graph_ids = sorted({
                    str(gene.get('graph_id', '') or '')
                    for _, gene in source_members
                    if str(gene.get('graph_id', '') or '')
                })
                reasons = sorted({
                    str(gene.get('reason', '') or '')
                    for member_status, gene in source_members
                    if member_status == 'unmatched' and gene.get('reason')
                })
                pathway_xrefs.append({
                    'database': database,
                    'gene_id': gene_id,
                    'gene_key': source_key,
                    'annotation_key': representative.get('annotation_key', ''),
                    'pathway_gene_name': names[0] if names else '',
                    'pathway_gene_names': names,
                    'graph_ids': graph_ids,
                    'node_occurrence_count': len(source_members),
                    'match_status': status,
                    'reason': ';'.join(reasons),
                    'match_method': representative.get('match_method', ''),
                    'mapping_source': representative.get('mapping_source', ''),
                    'matched_annotation_key': representative.get(
                        'matched_annotation_key', ''
                    ),
                    'matched_annotation_id': representative.get(
                        'matched_annotation_id', ''
                    ),
                })

            if matched_members:
                canonical_subject_ids = sorted({
                    str(gene.get('canonical_subject_id', '') or '')
                    for gene in matched_members
                    if str(gene.get('canonical_subject_id', '') or '')
                })
                if len(canonical_subject_ids) > 1:
                    raise ValueError(
                        "Conflicting canonical subjects for PW-KFI gene "
                        f"{pw_kfi_gene_key}: {canonical_subject_ids}"
                    )
                unit = min(
                    matched_members, key=self._pathway_member_sort_key
                ).copy()
                guides = sorted({
                    guide
                    for gene in matched_members
                    for guide in self._guide_entries(
                        gene.get('guide_targets', '')
                    )
                })
                query_ids = sorted({
                    query_id
                    for gene in matched_members
                    for query_id in gene.get('annotation_query_ids', [])
                    if query_id
                })
                annotation_subject_ids = sorted({
                    str(gene.get('resolved_annotation_subject_id', '') or '')
                    for gene in matched_members
                    if str(gene.get(
                        'resolved_annotation_subject_id', ''
                    ) or '')
                })
                annotation_unique_ids = sorted({
                    str(gene.get('unique_id', '') or '')
                    for gene in matched_members
                    if str(gene.get('unique_id', '') or '')
                })
                annotation_source_databases = sorted({
                    str(gene.get('annotation_source_database', '') or '')
                    for gene in matched_members
                    if str(gene.get(
                        'annotation_source_database', ''
                    ) or '')
                })
                annotation_identities = sorted({
                    (
                        str(gene.get(
                            'resolved_annotation_subject_id', ''
                        ) or ''),
                        str(gene.get('unique_id', '') or ''),
                        str(gene.get(
                            'annotation_source_database', ''
                        ) or ''),
                    )
                    for gene in matched_members
                })
                annotation_conflict = len(annotation_identities) > 1
                if annotation_conflict:
                    for xref in pathway_xrefs:
                        xref['match_status'] = 'unmatched'
                        xref['reason'] = (
                            'ambiguous_canonical_annotation_identity'
                        )
                    guides = []
                unit.update({
                    'guide_targets': ','.join(guides),
                    'guide_count': len(guides),
                    'has_guides': bool(guides),
                    'annotation_query_ids': query_ids,
                    'annotation_subject_ids': annotation_subject_ids,
                    'annotation_unique_ids': annotation_unique_ids,
                    'annotation_source_databases': (
                        annotation_source_databases
                    ),
                    'annotation_record_count': len(annotation_identities),
                    'annotation_identity_candidates': [
                        '|'.join(identity)
                        for identity in annotation_identities
                    ],
                    'reason': (
                        'ambiguous_canonical_annotation_identity'
                        if annotation_conflict else ''
                    ),
                    'match_status': (
                        'unmatched' if annotation_conflict else 'matched'
                    ),
                })
                destination = (
                    unmatched_units if annotation_conflict else matched_units
                )
            else:
                unit = members[0][1].copy()
                reasons = sorted({
                    str(gene.get('reason', '') or '')
                    for _, gene in members if gene.get('reason')
                })
                unit.update({
                    'reason': ';'.join(reasons),
                    'guide_targets': '',
                    'guide_count': 0,
                    'has_guides': False,
                    'match_status': 'unmatched',
                })
                destination = unmatched_units

            representative_xref = pathway_xrefs[0]
            pathway_gene_names = sorted({
                name
                for xref in pathway_xrefs
                for name in xref['pathway_gene_names']
            })
            canonicalization_methods = sorted({
                str(gene.get('canonicalization_method', '') or '')
                for _, gene in members
                if str(gene.get('canonicalization_method', '') or '')
            })
            unit.update({
                'pathway_gene_database': representative_xref['database'],
                'pathway_gene_id': representative_xref['gene_id'],
                'pathway_gene_name': (
                    pathway_gene_names[0] if pathway_gene_names else ''
                ),
                'pathway_gene_names': pathway_gene_names,
                'gene_key': representative_xref['gene_key'],
                'annotation_key': representative_xref['annotation_key'],
                'pw_kfi_gene_key': pw_kfi_gene_key,
                'pathway_xrefs': pathway_xrefs,
                'source_databases': sorted({
                    xref['database'] for xref in pathway_xrefs
                }),
                'node_occurrence_count': sum(
                    xref['node_occurrence_count'] for xref in pathway_xrefs
                ),
                'unique_xref_count': len(pathway_xrefs),
                'canonicalization_methods': canonicalization_methods,
                'canonicalization_complete': all(
                    bool(gene.get('canonicalization_complete', False))
                    for _, gene in members
                ),
            })
            destination.append(unit)

        return matched_units, unmatched_units
    
    def compute_pathway_kfi(self, pathway_info, matched_genes):
        """Calculate Pathway Knock-out Feasibility Index (PW-KFI)"""
        total_genes = pathway_info['total_genes']
        if total_genes == 0:
            return 0.0, "Not targetable", 0

        gene_keys = [gene.get('pw_kfi_gene_key', '') for gene in matched_genes]
        if any(not gene_key for gene_key in gene_keys):
            raise ValueError(
                'PW-KFI requires a canonical key for every matched gene'
            )
        if len(gene_keys) != len(set(gene_keys)):
            raise ValueError(
                'PW-KFI matched genes must already be canonicalized and unique'
            )
        if len(matched_genes) > total_genes:
            raise ValueError(
                'Matched PW-KFI genes cannot exceed the pathway denominator'
            )
        
        genes_with_guides = sum(1 for gene in matched_genes if gene['has_guides'])
        pw_kfi_score = genes_with_guides / total_genes
        
        # 分類判定
        if pw_kfi_score == 1.0:
            targetable_type = "Fully targetable"
        elif pw_kfi_score > 0:
            targetable_type = "Partially targetable"  
        else:
            targetable_type = "Not targetable"
        
        return pw_kfi_score, targetable_type, genes_with_guides

    @staticmethod
    def _pathway_gene_audit_record(
        pathway_info, gene, status, reason='', included_in_denominator=None
    ):
        """Build one provenance-preserving pathway gene audit row."""
        database = gene.get('pathway_gene_database', gene.get('database', ''))
        gene_id = gene.get('pathway_gene_id', gene.get('gene_id', ''))
        if included_in_denominator is None:
            included_in_denominator = status != 'not_evaluable'
        return {
            'pathway_id': pathway_info.get('pathway_id', ''),
            'pathway_name': pathway_info.get('pathway_name', ''),
            'gpml_file_path': pathway_info.get('gpml_file_path', ''),
            'pathway_gene_database': database,
            'pathway_gene_id': gene_id,
            'pathway_gene_name': gene.get(
                'pathway_gene_name', gene.get('gene_name', '')
            ),
            'gene_key': f"{database}::{gene_id}" if database or gene_id else '',
            'annotation_key': gene.get('annotation_key', ''),
            'pw_kfi_gene_key': gene.get('pw_kfi_gene_key', ''),
            'node_occurrence_count': int(gene.get(
                'node_occurrence_count', 1
            )),
            'denominator_contribution': int(bool(included_in_denominator)),
            'included_in_denominator': bool(included_in_denominator),
            'status': status,
            'reason': reason,
            'canonical_subject_id': gene.get('canonical_subject_id', ''),
            'canonicalization_methods': ','.join(
                gene.get('canonicalization_methods', [])
            ) if isinstance(
                gene.get('canonicalization_methods', []), list
            ) else gene.get('canonicalization_methods', ''),
            'canonicalization_complete': bool(gene.get(
                'canonicalization_complete', False
            )),
            'pathway_canonicalization_status': pathway_info.get(
                'canonicalization_status', ''
            ),
            'match_method': gene.get('match_method', ''),
            'mapping_source': gene.get('mapping_source', ''),
            'annotation_source_database': gene.get(
                'annotation_source_database', ''
            ),
            'matched_annotation_key': gene.get('matched_annotation_key', ''),
            'matched_annotation_id': gene.get('matched_annotation_id', ''),
            'resolved_annotation_subject_id': gene.get(
                'resolved_annotation_subject_id', ''
            ),
            'sgd_id': gene.get('sgd_id', ''),
            'ensembl_id': gene.get('ensembl_id', ''),
            'unique_id': gene.get('unique_id', ''),
        }

    def _audit_collapsed_gene_unit(self, pathway_info, gene_unit):
        """Audit every exact source xref while counting one canonical gene."""
        for index, xref in enumerate(gene_unit.get('pathway_xrefs', [])):
            audit_gene = gene_unit.copy()
            audit_gene.update({
                'pathway_gene_database': xref['database'],
                'pathway_gene_id': xref['gene_id'],
                'pathway_gene_name': xref.get('pathway_gene_name', ''),
                'gene_key': xref['gene_key'],
                'annotation_key': xref.get('annotation_key', ''),
                'node_occurrence_count': xref.get('node_occurrence_count', 1),
                'match_method': xref.get('match_method', ''),
                'mapping_source': xref.get('mapping_source', ''),
                'matched_annotation_key': xref.get(
                    'matched_annotation_key', ''
                ),
                'matched_annotation_id': xref.get(
                    'matched_annotation_id', ''
                ),
            })
            status = xref.get('match_status', gene_unit['match_status'])
            reason = xref.get('reason', '') or (
                gene_unit.get('reason', '') if status == 'unmatched' else ''
            )
            self.pathway_gene_audit.append(
                self._pathway_gene_audit_record(
                    pathway_info,
                    audit_gene,
                    status=status,
                    reason=reason,
                    included_in_denominator=index == 0,
                )
            )
    
    def process_all_pathways(self):
        """Process all pathways"""
        self.logger.info("Analysis of all pathways has begun...")

        self.pathway_gene_audit = []
        self.coverage_stats = {
            'total_pathways': 0,
            'processable_pathways': 0,
            'skipped_pathways': 0,
            'coverage_by_database': {},
        }
        
        gpml_files = self.get_gpml_files()
        annotation_data = self.load_annotation_results()
        
        processed_pathways = []
        
        for gpml_file in gpml_files:
            self.logger.info(f"Processing: {gpml_file.name}")
            
            # GPMLファイル解析
            pathway_info = self.parse_gpml_file(gpml_file)
            if not pathway_info:
                continue

            for rejected_gene in pathway_info.get('unmatchable_genes', []):
                self.pathway_gene_audit.append(
                    self._pathway_gene_audit_record(
                        pathway_info,
                        rejected_gene,
                        status='not_evaluable',
                        reason=rejected_gene.get('reason', 'not_evaluable'),
                        included_in_denominator=False,
                    )
                )
                database = rejected_gene.get('database', '') or 'Unspecified'
                database_stats = self.coverage_stats['coverage_by_database'].setdefault(
                    database,
                    {'evaluable': 0, 'matched': 0, 'unmatched': 0, 'not_evaluable': 0},
                )
                database_stats['not_evaluable'] += 1

            if not pathway_info['genes']:
                continue

            unique_source_xrefs = {
                (gene['database'], gene['gene_id'])
                for gene in pathway_info['genes']
            }
            for database, _ in sorted(unique_source_xrefs):
                database_stats = self.coverage_stats['coverage_by_database'].setdefault(
                    database,
                    {'evaluable': 0, 'matched': 0, 'unmatched': 0, 'not_evaluable': 0},
                )
                database_stats['evaluable'] += 1
            
            # 遺伝子照合
            matched_genes, unmatched_genes = self.match_pathway_genes_with_annotation(
                pathway_info['genes'], annotation_data
            )
            matched_genes, unmatched_genes = self.collapse_pathway_gene_units(
                matched_genes, unmatched_genes
            )

            pathway_info['total_genes'] = (
                len(matched_genes) + len(unmatched_genes)
            )
            pathway_info['gene_count_basis'] = 'exact_canonical_gene_units'
            source_databases = {
                gene['database'] for gene in pathway_info['genes']
            }
            all_gene_units = matched_genes + unmatched_genes
            canonicalization_complete = (
                len(source_databases) <= 1
                or all(
                    gene.get('canonicalization_complete', False)
                    for gene in all_gene_units
                )
            )
            pathway_info['canonicalization_complete'] = (
                canonicalization_complete
            )
            pathway_info['canonicalization_status'] = (
                'complete'
                if canonicalization_complete
                else 'incomplete_cross_database_identity'
            )

            for gene_unit in all_gene_units:
                self._audit_collapsed_gene_unit(pathway_info, gene_unit)
                for xref in gene_unit['pathway_xrefs']:
                    database_stats = self.coverage_stats[
                        'coverage_by_database'
                    ][xref['database']]
                    database_stats[xref['match_status']] += 1

            if not canonicalization_complete:
                self.logger.warning(
                    "PW-KFI calculation skipped because exact cross-database "
                    "gene identity could not be resolved: "
                    f"{pathway_info['pathway_name']}"
                )
                self.coverage_stats['skipped_pathways'] += 1
                continue
            
            # カバー率計算
            coverage_rate = len(matched_genes) / pathway_info['total_genes']
            min_coverage = self.pathway_config.get('minimum_pathway_coverage', 0.5)
            
            if coverage_rate < min_coverage:
                self.logger.warning(
                    f"Low coverage pathways: {pathway_info['pathway_name']} "
                    f"({coverage_rate:.1%} < {min_coverage:.1%})"
                )
                if self.pathway_config.get('skip_low_coverage_pathways', True):
                    self.coverage_stats['skipped_pathways'] += 1
                    continue

            if not matched_genes:
                self.logger.warning(
                    f"No exact annotation match: {pathway_info['pathway_name']}"
                )
            
            # PW-KFI計算
            pw_kfi_score, targetable_type, genes_with_guides = self.compute_pathway_kfi(
                pathway_info, matched_genes
            )
            
            # 結果統合
            pathway_result = {
                **pathway_info,
                'matched_genes': matched_genes,
                'unmatched_genes': unmatched_genes,
                'matched_gene_count': len(matched_genes),
                'unmatched_gene_count': len(unmatched_genes),
                'genes_with_guides': genes_with_guides,
                'pw_kfi_score': pw_kfi_score,
                'targetable_type': targetable_type,
                'coverage_rate': coverage_rate,
                'match_counts_by_database': dict(Counter(
                    xref['database']
                    for gene in all_gene_units
                    for xref in gene['pathway_xrefs']
                    if xref['match_status'] == 'matched'
                )),
                'unmatched_counts_by_database': dict(Counter(
                    xref['database']
                    for gene in all_gene_units
                    for xref in gene['pathway_xrefs']
                    if xref['match_status'] == 'unmatched'
                )),
            }
            
            processed_pathways.append(pathway_result)
            
            self.logger.info(
                f"完了: {pathway_info['pathway_name']} - "
                f"PW-KFI: {pw_kfi_score:.3f} ({targetable_type})"
            )
        
        return processed_pathways
    
    def save_pathway_results(self, processed_pathways):
        """Save pathway analysis results"""
        # PW-KFIテーブル保存
        kfi_table_path = os.path.join(self.intermediate_dir, 'pathway_kfi_table.csv')
        
        with open(kfi_table_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = [
                'pathway_id', 'pathway_name', 'pathway_organism',
                'total_gene_nodes', 'evaluable_gene_node_count',
                'unique_xref_count', 'total_genes', 'gene_count_basis',
                'canonicalization_complete', 'canonicalization_status',
                'unmatchable_gene_count',
                'unmatchable_genes', 'matched_genes', 'unmatched_genes',
                'matched_gene_count',
                'unmatched_gene_count', 'genes_with_guides', 'pw_kfi_score',
                'targetable_type', 'coverage_rate', 'match_counts_by_database',
                'unmatched_counts_by_database', 'gpml_file_path'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_NONNUMERIC)
            writer.writeheader()
            
            for pathway in processed_pathways:
                row = {key: pathway.get(key, '') for key in fieldnames}
                writer.writerow(row)
        
        self.logger.info(f"PW-KFI Table Storage: {kfi_table_path}")
        
        # 完全標的可能パスウェイ抽出
        fully_targetable = [p for p in processed_pathways if p['targetable_type'] == 'Fully targetable']
        fully_targetable_path = os.path.join(self.intermediate_dir, 'pathway_kfi_fully_targetable.csv')
        
        with open(fully_targetable_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = [
                'pathway_id', 'pathway_name', 'total_genes', 'pw_kfi_score', 'coverage_rate'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames, quoting=csv.QUOTE_NONNUMERIC)
            writer.writeheader()
            
            for pathway in fully_targetable:
                row = {key: pathway.get(key, '') for key in fieldnames}
                writer.writerow(row)
        
        self.logger.info(f"Fully targetable pathway preservation: {fully_targetable_path} ({len(fully_targetable)}個)")
        
        # パスウェイ-遺伝子マッピング保存
        mapping_path = os.path.join(self.intermediate_dir, 'pathway_gene_mapping.tsv')
        
        with open(mapping_path, 'w', encoding='utf-8', newline='') as f:
            fieldnames = [
                'pathway_id', 'pathway_name', 'pathway_gene_database',
                'pathway_gene_id', 'pathway_gene_name', 'gene_key',
                'pw_kfi_gene_key', 'pathway_xrefs', 'source_databases',
                'node_occurrence_count', 'unique_xref_count',
                'canonical_subject_id', 'canonicalization_methods',
                'canonicalization_complete',
                'annotation_source_database',
                'matched_annotation_key', 'matched_annotation_id',
                'resolved_annotation_subject_id', 'match_method', 'mapping_source',
                'sgd_id', 'ensembl_id', 'gene_sgd_id', 'unique_id', 'gene_name',
                'has_guides', 'guide_count',
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', quoting=csv.QUOTE_NONNUMERIC)
            writer.writeheader()
            
            for pathway in processed_pathways:
                for gene in pathway['matched_genes']:
                    # gene_name の安全な取得
                    gene_name = gene.get('gene_name', 'Unknown')
                    if not isinstance(gene_name, str):
                        gene_name = str(gene_name)
                    
                    writer.writerow({
                        'pathway_id': str(pathway['pathway_id']),
                        'pathway_name': str(pathway['pathway_name']),
                        'pathway_gene_database': gene.get('pathway_gene_database', ''),
                        'pathway_gene_id': gene.get('pathway_gene_id', ''),
                        'pathway_gene_name': gene.get('pathway_gene_name', ''),
                        'gene_key': gene.get('gene_key', ''),
                        'pw_kfi_gene_key': gene.get('pw_kfi_gene_key', ''),
                        'pathway_xrefs': repr(gene.get('pathway_xrefs', [])),
                        'source_databases': ','.join(
                            gene.get('source_databases', [])
                        ),
                        'node_occurrence_count': gene.get(
                            'node_occurrence_count', 0
                        ),
                        'unique_xref_count': gene.get('unique_xref_count', 0),
                        'canonical_subject_id': gene.get(
                            'canonical_subject_id', ''
                        ),
                        'canonicalization_methods': ','.join(
                            gene.get('canonicalization_methods', [])
                        ),
                        'canonicalization_complete': bool(gene.get(
                            'canonicalization_complete', False
                        )),
                        'annotation_source_database': gene.get(
                            'annotation_source_database', ''
                        ),
                        'matched_annotation_key': gene.get('matched_annotation_key', ''),
                        'matched_annotation_id': gene.get('matched_annotation_id', ''),
                        'resolved_annotation_subject_id': gene.get(
                            'resolved_annotation_subject_id', ''
                        ),
                        'match_method': gene.get('match_method', ''),
                        'mapping_source': gene.get('mapping_source', ''),
                        'sgd_id': gene.get('sgd_id', ''),
                        'ensembl_id': gene.get('ensembl_id', ''),
                        'gene_sgd_id': gene.get('sgd_id', ''),
                        'unique_id': gene.get('unique_id', ''),
                        'gene_name': gene_name,
                        'has_guides': bool(gene['has_guides']),
                        'guide_count': int(gene['guide_count'])
                    })
        
        self.logger.info(f"Pathway-gene mapping preservation: {mapping_path}")

        audit_path = os.path.join(
            self.intermediate_dir, 'pathway_gene_match_audit.tsv'
        )
        audit_fieldnames = [
            'pathway_id', 'pathway_name', 'gpml_file_path',
            'pathway_gene_database', 'pathway_gene_id', 'pathway_gene_name',
            'gene_key', 'pw_kfi_gene_key', 'annotation_key',
            'node_occurrence_count', 'denominator_contribution',
            'included_in_denominator', 'status',
            'reason', 'match_method', 'mapping_source',
            'canonical_subject_id', 'canonicalization_methods',
            'canonicalization_complete', 'pathway_canonicalization_status',
            'annotation_source_database', 'matched_annotation_key',
            'matched_annotation_id', 'resolved_annotation_subject_id',
            'sgd_id', 'ensembl_id', 'unique_id',
        ]
        with open(audit_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(
                f, fieldnames=audit_fieldnames, delimiter='\t',
                quoting=csv.QUOTE_NONNUMERIC,
            )
            writer.writeheader()
            writer.writerows(self.pathway_gene_audit)
        self.logger.info(f"Pathway gene match audit preservation: {audit_path}")

        annotation_audit_path = os.path.join(
            self.intermediate_dir, 'pathway_annotation_identifier_audit.tsv'
        )
        annotation_audit_fieldnames = [
            'source_path', 'row_number', 'database', 'annotation_key',
            'identifier', 'query_id', 'unique_id', 'reason',
            'conflicting_fields', 'candidate_values',
        ]
        with open(annotation_audit_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(
                f, fieldnames=annotation_audit_fieldnames, delimiter='\t',
                quoting=csv.QUOTE_NONNUMERIC, extrasaction='ignore',
            )
            writer.writeheader()
            writer.writerows(self.annotation_identifier_issues)
        self.logger.info(
            f"Annotation identifier audit preservation: {annotation_audit_path}"
        )
        
        return {
            'kfi_table_path': kfi_table_path,
            'fully_targetable_path': fully_targetable_path,
            'mapping_path': mapping_path,
            'gene_match_audit_path': audit_path,
            'annotation_identifier_audit_path': annotation_audit_path,
        }
    
    def generate_statistics_report(self, processed_pathways):
        """Generate statistics report"""
        audit_status_counts = Counter(
            row['status'] for row in self.pathway_gene_audit
        )
        audit_reason_counts = Counter(
            row['reason'] for row in self.pathway_gene_audit if row['reason']
        )
        stats = {
            'total_pathways': self.coverage_stats['total_pathways'],
            'processable_pathways': len(processed_pathways),
            'skipped_pathways': self.coverage_stats['skipped_pathways'],
            'fully_targetable_count': len([p for p in processed_pathways if p['targetable_type'] == 'Fully targetable']),
            'partially_targetable_count': len([p for p in processed_pathways if p['targetable_type'] == 'Partially targetable']),
            'not_targetable_count': len([p for p in processed_pathways if p['targetable_type'] == 'Not targetable']),
            'pathway_gene_status_counts': dict(audit_status_counts),
            'pathway_gene_reason_counts': dict(audit_reason_counts),
            'coverage_by_database': self.coverage_stats['coverage_by_database'],
            'annotation_identifier_issue_count': len(self.annotation_identifier_issues),
        }
        
        # パーセンテージ計算
        if stats['processable_pathways'] > 0:
            stats['fully_targetable_percent'] = (stats['fully_targetable_count'] / stats['processable_pathways']) * 100
            stats['partially_targetable_percent'] = (stats['partially_targetable_count'] / stats['processable_pathways']) * 100
            stats['not_targetable_percent'] = (stats['not_targetable_count'] / stats['processable_pathways']) * 100
        
        # 統計レポート保存
        stats_path = os.path.join(self.intermediate_dir, 'pathway_analysis_results.txt')
        
        with open(stats_path, 'w', encoding='utf-8') as f:
            f.write("=== Pathway Analysis Statistical Report ===\n\n")
            f.write(f"Total number of pathways: {stats['total_pathways']}\n")
            total_pathways = stats['total_pathways']
            processable_percent = (
                stats['processable_pathways'] / total_pathways * 100
                if total_pathways else 0.0
            )
            skipped_percent = (
                stats['skipped_pathways'] / total_pathways * 100
                if total_pathways else 0.0
            )
            f.write(f"Number of processable pathways: {stats['processable_pathways']} ({processable_percent:.1f}%)\n")
            f.write(f"Number of skipped pathways: {stats['skipped_pathways']} ({skipped_percent:.1f}%)\n\n")
            
            if stats['processable_pathways'] > 0:
                f.write("=== PW-KFI classification statistics ===\n")
                f.write(f"Fully targetable: {stats['fully_targetable_count']} ({stats['fully_targetable_percent']:.1f}%)\n")
                f.write(f"Partial targeting possible: {stats['partially_targetable_count']} ({stats['partially_targetable_percent']:.1f}%)\n")
                f.write(f"Untargetable: {stats['not_targetable_count']} ({stats['not_targetable_percent']:.1f}%)\n\n")

            f.write("=== Exact identifier matching by source database ===\n")
            for database in sorted(stats['coverage_by_database']):
                database_stats = stats['coverage_by_database'][database]
                f.write(
                    f"{database}: evaluable={database_stats['evaluable']}, "
                    f"matched={database_stats['matched']}, "
                    f"unmatched={database_stats['unmatched']}, "
                    f"not_evaluable={database_stats['not_evaluable']}\n"
                )

            f.write("\n=== Unmatched and non-evaluable reasons ===\n")
            if audit_reason_counts:
                for reason in sorted(audit_reason_counts):
                    f.write(f"{reason}: {audit_reason_counts[reason]}\n")
            else:
                f.write("None\n")
            f.write(
                "Annotation identifier audit issues: "
                f"{len(self.annotation_identifier_issues)}\n\n"
            )
            
            f.write("=== Configuration Information ===\n")
            f.write(f"GPML directory: {self.pathway_config['wikipathways_gpml_dir']}\n")
            f.write(f"Target database: {[db['database_name'] for db in self.pathway_config['target_gene_databases']]}\n")
            f.write(f"Minimum coverage: {self.pathway_config.get('minimum_pathway_coverage', 0.5)}\n")
        
        self.logger.info(f"Statistical report storage: {stats_path}")
        return stats
    
    def run(self):
        """Execute entire pathway analysis"""
        try:
            self.logger.info("Start pathway analysis...")
            
            # 設定検証
            self.validate_configuration()
            
            # 全パスウェイ処理
            processed_pathways = self.process_all_pathways()

            if not processed_pathways:
                self.logger.warning(
                    "No processable pathways were found; audit artifacts will still be saved."
                )
            
            # 結果保存
            saved_files = self.save_pathway_results(processed_pathways)
            
            # Generate statistics report
            stats = self.generate_statistics_report(processed_pathways)
            
            # 結果まとめ
            results = {
                'processed_pathways': processed_pathways,
                'statistics': stats,
                'output_files': saved_files
            }
            
            self.logger.info(f"Pathway analysis completed: {len(processed_pathways)} processed")
            return results
            
        except Exception as e:
            self.logger.error(f"An error occurred during pathway analysis.: {e}")
            raise
