import csv
import tempfile
import unittest
from pathlib import Path

from src.modules.pathway_analysis import PathwayAnalysis
from src.modules.pathway_report_generator import PathwayReportGenerator


REPO_ROOT = Path(__file__).resolve().parents[1]
BUNDLED_GPML_DIR = (
    REPO_ROOT
    / 'LOCAL_DATA/db/wikipathways-20250910-gpml-Saccharomyces_cerevisiae'
)
BUNDLED_SGD_PATH = (
    REPO_ROOT / 'LOCAL_DATA/db/sgd_S288C_reference_orf_trans_headermodified.fasta'
)


ANNOTATION_COLUMNS = [
    'query_id', 'subject_id', 'identity', 'length', 'evalue', 'bit_score',
    'description', 'gene_name', 'classification', 'unique_id',
    'annotation_source_database', 'sgd_id', 'ensembl_id', 'guide_targets',
]


class PathwayAnalysisIdentifierTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir_context = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temp_dir_context.name)
        self.gpml_dir = self.temp_dir / 'gpml'
        self.gpml_dir.mkdir()
        self.sgd_path = self.temp_dir / 'sgd.fasta'
        self.sgd_path.write_text(
            '>YAL060W|BDH1|SGDID:S000000056 NAD-dependent protein\nMPEPTIDE\n',
            encoding='utf-8',
        )
        self.config = {
            'input': {'sgd_db': str(self.sgd_path)},
            'output': {
                'intermediate': str(self.temp_dir),
                'base_dir': str(self.temp_dir / 'output'),
            },
            'pathway_analysis': {
                'wikipathways_gpml_dir': str(self.gpml_dir),
                'target_gene_databases': [
                    {
                        'database_name': 'SGD',
                        'id_pattern': r'^S\d{9}$',
                        'annotation_key': 'sgd_id',
                    },
                    {
                        'database_name': 'Ensembl',
                        'id_pattern': r'^Y[A-Z]{2}\d{3}[CW]$',
                        'annotation_key': 'ensembl_id',
                    },
                ],
                'minimum_pathway_coverage': 0.5,
                'skip_low_coverage_pathways': False,
            },
        }

    def tearDown(self):
        self.temp_dir_context.cleanup()

    def write_annotation(self, rows, columns=ANNOTATION_COLUMNS):
        path = self.temp_dir / 'annotation_with_guides_final.tsv'
        with path.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns, delimiter='\t', extrasaction='ignore'
            )
            writer.writeheader()
            writer.writerows(rows)
        return path

    def write_gpml(self, nodes, name='Strict pathway', pathway_id='WP9999'):
        node_xml = []
        for index, (database, gene_id, label) in enumerate(nodes, 1):
            node_xml.append(
                f'<DataNode TextLabel="{label}" GraphId="n{index}" '
                f'Type="GeneProduct"><Xref Database="{database}" '
                f'ID="{gene_id}" /></DataNode>'
            )
        path = self.gpml_dir / f'Sc_test_{pathway_id}_20260715.gpml'
        path.write_text(
            '<?xml version="1.0" encoding="UTF-8"?>\n'
            '<Pathway xmlns="http://pathvisio.org/GPML/2013a" '
            f'Name="{name}" Organism="Saccharomyces cerevisiae" '
            f'Version="20260715 {pathway_id}">'
            + ''.join(node_xml)
            + '</Pathway>',
            encoding='utf-8',
        )
        return path

    @staticmethod
    def annotation_row(**overrides):
        row = {
            'query_id': 'query-1',
            'subject_id': 'SGDID:S000000056',
            'identity': '',
            'length': '383',
            'evalue': '0',
            'bit_score': '154.8',
            'description': 'NAD-dependent protein',
            'gene_name': 'BDH1',
            'classification': 'Single',
            'unique_id': 'Single_0000000001',
            'annotation_source_database': 'SGD',
            'sgd_id': 'S000000056',
            'ensembl_id': 'YAL060W',
            'guide_targets': 'g1/1/23/ACGT/0/0',
        }
        row.update(overrides)
        return row

    def test_matches_sgd_and_ensembl_by_separate_exact_identifiers(self):
        self.write_annotation([self.annotation_row()])
        analysis = PathwayAnalysis(self.config)
        annotation_data = analysis.load_annotation_results()

        matched, unmatched = analysis.match_pathway_genes_with_annotation(
            [
                {
                    'database': 'SGD',
                    'gene_id': 'S000000056',
                    'gene_name': 'BDH1',
                    'annotation_key': 'sgd_id',
                },
                {
                    'database': 'Ensembl',
                    'gene_id': 'YAL060W',
                    'gene_name': 'BDH1',
                    'annotation_key': 'ensembl_id',
                },
            ],
            annotation_data,
        )

        self.assertEqual(unmatched, [])
        self.assertEqual(len(matched), 2)
        self.assertEqual(
            {row['pathway_gene_database'] for row in matched},
            {'SGD', 'Ensembl'},
        )
        ensembl_match = next(
            row for row in matched if row['pathway_gene_database'] == 'Ensembl'
        )
        self.assertEqual(ensembl_match['pathway_gene_id'], 'YAL060W')
        self.assertEqual(ensembl_match['annotation_source_database'], 'SGD')
        self.assertEqual(ensembl_match['matched_annotation_key'], 'ensembl_id')
        self.assertEqual(ensembl_match['matched_annotation_id'], 'YAL060W')
        self.assertEqual(ensembl_match['match_method'], 'exact_annotation_identifier')

    def test_uses_exact_sgd_fasta_cross_reference_but_not_gene_name(self):
        legacy_columns = [
            'query_id', 'subject_id', 'identity', 'length', 'evalue',
            'bit_score', 'description', 'gene_name', 'classification',
            'unique_id', 'guide_targets',
        ]
        self.write_annotation(
            [self.annotation_row(gene_name='YBR160W')],
            columns=legacy_columns,
        )
        analysis = PathwayAnalysis(self.config)
        annotation_data = analysis.load_annotation_results()

        matched, unmatched = analysis.match_pathway_genes_with_annotation(
            [
                {
                    'database': 'Ensembl',
                    'gene_id': 'YAL060W',
                    'gene_name': 'anything',
                    'annotation_key': 'ensembl_id',
                },
                {
                    'database': 'Ensembl',
                    'gene_id': 'YBR160W',
                    'gene_name': 'YBR160W',
                    'annotation_key': 'ensembl_id',
                },
            ],
            annotation_data,
        )

        self.assertEqual(len(matched), 1)
        self.assertEqual(matched[0]['pathway_gene_id'], 'YAL060W')
        self.assertEqual(
            matched[0]['match_method'], 'exact_sgd_fasta_cross_reference'
        )
        self.assertEqual(matched[0]['annotation_source_database'], 'SGD')
        self.assertEqual(len(unmatched), 1)
        self.assertEqual(unmatched[0]['gene_id'], 'YBR160W')
        self.assertEqual(unmatched[0]['reason'], 'cross_reference_not_found')

    def test_rejects_ambiguous_annotation_identifier(self):
        self.write_annotation([
            self.annotation_row(),
            self.annotation_row(
                query_id='query-2',
                subject_id='SGDID:S000000057',
                sgd_id='S000000057',
                unique_id='Single_0000000002',
            ),
        ])
        analysis = PathwayAnalysis(self.config)
        annotation_data = analysis.load_annotation_results()

        matched, unmatched = analysis.match_pathway_genes_with_annotation(
            [{
                'database': 'Ensembl',
                'gene_id': 'YAL060W',
                'gene_name': 'BDH1',
                'annotation_key': 'ensembl_id',
            }],
            annotation_data,
        )

        self.assertEqual(matched, [])
        self.assertEqual(unmatched[0]['reason'], 'ambiguous_annotation_identifier')

    def test_rejects_ambiguous_exact_sgd_fasta_cross_reference(self):
        self.sgd_path.write_text(
            '>YAL060W|BDH1|SGDID:S000000056 first\nMPEPTIDE\n'
            '>YAL060W|BDH1|SGDID:S000000057 conflicting\nMPEPTIDE\n',
            encoding='utf-8',
        )
        legacy_columns = [
            'query_id', 'subject_id', 'identity', 'length', 'evalue',
            'bit_score', 'description', 'gene_name', 'classification',
            'unique_id', 'guide_targets',
        ]
        self.write_annotation([self.annotation_row()], columns=legacy_columns)
        analysis = PathwayAnalysis(self.config)
        annotation_data = analysis.load_annotation_results()

        matched, unmatched = analysis.match_pathway_genes_with_annotation(
            [{
                'database': 'Ensembl',
                'gene_id': 'YAL060W',
                'gene_name': 'BDH1',
                'annotation_key': 'ensembl_id',
            }],
            annotation_data,
        )

        self.assertEqual(matched, [])
        self.assertEqual(unmatched[0]['reason'], 'ambiguous_cross_reference')

    def test_rejects_sgd_annotation_alias_conflicting_with_cross_reference(self):
        self.write_annotation([
            self.annotation_row(
                subject_id='SGDID:S000000057',
                sgd_id='S000000057',
            )
        ])
        analysis = PathwayAnalysis(self.config)
        annotation_data = analysis.load_annotation_results()

        matched, unmatched = analysis.match_pathway_genes_with_annotation(
            [{
                'database': 'Ensembl',
                'gene_id': 'YAL060W',
                'gene_name': 'BDH1',
                'annotation_key': 'ensembl_id',
            }],
            annotation_data,
        )

        self.assertEqual(matched, [])
        self.assertEqual(unmatched[0]['reason'], 'conflicting_cross_reference')

    def test_rejects_subject_id_conflicting_with_explicit_sgd_id(self):
        self.write_annotation([
            self.annotation_row(
                subject_id='SGDID:S000000056',
                sgd_id='S000000057',
                ensembl_id='',
            )
        ])
        analysis = PathwayAnalysis(self.config)
        annotation_data = analysis.load_annotation_results()

        matched, unmatched = analysis.match_pathway_genes_with_annotation(
            [{
                'database': 'SGD',
                'gene_id': 'S000000056',
                'gene_name': 'BDH1',
                'annotation_key': 'sgd_id',
            }],
            annotation_data,
        )

        self.assertEqual(matched, [])
        self.assertEqual(
            unmatched[0]['reason'], 'ambiguous_annotation_subject_id'
        )
        self.assertIn(
            'SGDID:S000000056', annotation_data['ambiguous_subject_ids']
        )
        self.assertEqual(
            [issue['reason'] for issue in analysis.annotation_identifier_issues],
            ['conflicting_annotation_identity'],
        )

    def test_direct_explicit_identifier_does_not_require_optional_resolver(self):
        self.write_annotation([self.annotation_row()])
        self.config['input']['sgd_db'] = str(self.temp_dir / 'missing.fasta')
        analysis = PathwayAnalysis(self.config)
        annotation_data = analysis.load_annotation_results()

        matched, unmatched = analysis.match_pathway_genes_with_annotation(
            [{
                'database': 'Ensembl',
                'gene_id': 'YAL060W',
                'gene_name': 'BDH1',
                'annotation_key': 'ensembl_id',
            }],
            annotation_data,
        )

        self.assertEqual(len(matched), 1)
        self.assertEqual(unmatched, [])
        self.assertEqual(matched[0]['match_method'], 'exact_annotation_identifier')

    def test_legacy_annotation_reports_unavailable_resolver_without_name_fallback(self):
        legacy_columns = [
            'query_id', 'subject_id', 'identity', 'length', 'evalue',
            'bit_score', 'description', 'gene_name', 'classification',
            'unique_id', 'guide_targets',
        ]
        self.write_annotation(
            [self.annotation_row(gene_name='YAL060W')],
            columns=legacy_columns,
        )
        self.config['input'].pop('sgd_db')
        analysis = PathwayAnalysis(self.config)
        annotation_data = analysis.load_annotation_results()

        matched, unmatched = analysis.match_pathway_genes_with_annotation(
            [{
                'database': 'Ensembl',
                'gene_id': 'YAL060W',
                'gene_name': 'YAL060W',
                'annotation_key': 'ensembl_id',
            }],
            annotation_data,
        )

        self.assertEqual(matched, [])
        self.assertEqual(
            unmatched[0]['reason'], 'resolver_source_not_configured'
        )

    def test_extracts_only_fullmatch_ids_and_reports_rejections(self):
        gpml_path = self.write_gpml([
            ('SGD', 'S000000056', 'valid SGD'),
            ('Ensembl', 'YAL060W', 'valid Ensembl'),
            ('SGD', 'BDH1', 'invalid SGD'),
            ('Ensembl', 'YAL060W-extra', 'invalid Ensembl'),
            ('Entrez Gene', '123', 'unsupported'),
        ])
        analysis = PathwayAnalysis(self.config)

        pathway = analysis.parse_gpml_file(gpml_path)

        self.assertEqual(pathway['total_genes'], 2)
        self.assertEqual(pathway['unmatchable_gene_count'], 3)
        self.assertEqual(
            [row['reason'] for row in pathway['unmatchable_genes']],
            [
                'invalid_identifier_format',
                'invalid_identifier_format',
                'unsupported_database',
            ],
        )

    def test_collapses_duplicate_nodes_and_cross_database_aliases_for_pw_kfi(self):
        self.write_annotation([self.annotation_row()])
        self.write_gpml([
            ('Ensembl', 'YAL060W', 'BDH1 first node'),
            ('Ensembl', 'YAL060W', 'BDH1 duplicate node'),
            ('SGD', 'S000000056', 'BDH1 SGD alias'),
        ])
        analysis = PathwayAnalysis(self.config)

        pathways = analysis.process_all_pathways()

        self.assertEqual(len(pathways), 1)
        pathway = pathways[0]
        self.assertEqual(pathway['total_gene_nodes'], 3)
        self.assertEqual(pathway['evaluable_gene_node_count'], 3)
        self.assertEqual(pathway['unique_xref_count'], 2)
        self.assertEqual(pathway['total_genes'], 1)
        self.assertEqual(pathway['matched_gene_count'], 1)
        self.assertEqual(pathway['unmatched_gene_count'], 0)
        self.assertEqual(pathway['genes_with_guides'], 1)
        self.assertEqual(pathway['pw_kfi_score'], 1.0)
        gene = pathway['matched_genes'][0]
        self.assertEqual(gene['pw_kfi_gene_key'], 'SGD::SGDID:S000000056')
        self.assertEqual(gene['node_occurrence_count'], 3)
        self.assertEqual(gene['unique_xref_count'], 2)
        self.assertEqual(
            [xref['gene_key'] for xref in gene['pathway_xrefs']],
            ['Ensembl::YAL060W', 'SGD::S000000056'],
        )
        self.assertEqual(
            [xref['node_occurrence_count'] for xref in gene['pathway_xrefs']],
            [2, 1],
        )
        self.assertEqual(
            sum(row['denominator_contribution'] for row in analysis.pathway_gene_audit),
            1,
        )

    def test_custom_annotation_subject_does_not_replace_exact_pathway_identity(self):
        self.write_annotation([
            self.annotation_row(
                subject_id='CUSTOM:BDH1',
                annotation_source_database='Custom',
                sgd_id='',
            )
        ])
        self.write_gpml([
            ('Ensembl', 'YAL060W', 'BDH1 Ensembl'),
            ('SGD', 'S000000056', 'BDH1 SGD'),
        ])
        analysis = PathwayAnalysis(self.config)

        pathway = analysis.process_all_pathways()[0]

        self.assertEqual(pathway['total_genes'], 1)
        self.assertEqual(pathway['matched_gene_count'], 1)
        self.assertEqual(pathway['genes_with_guides'], 1)
        self.assertEqual(pathway['pw_kfi_score'], 1.0)
        gene = pathway['matched_genes'][0]
        self.assertEqual(gene['canonical_subject_id'], 'SGDID:S000000056')
        self.assertEqual(
            gene['resolved_annotation_subject_id'], 'CUSTOM:BDH1'
        )
        self.assertEqual(
            gene['annotation_subject_ids'], ['CUSTOM:BDH1']
        )
        self.assertEqual(
            [xref['match_status'] for xref in gene['pathway_xrefs']],
            ['matched', 'unmatched'],
        )

    def test_conflicting_exact_alias_annotations_fail_closed_per_gene(self):
        self.write_annotation([
            self.annotation_row(
                query_id='custom-query',
                subject_id='CUSTOM:BDH1',
                unique_id='Custom_0001',
                annotation_source_database='Custom',
                sgd_id='',
                ensembl_id='YAL060W',
                guide_targets='custom/1/23/ACGT/0/0',
            ),
            self.annotation_row(
                query_id='sgd-query',
                subject_id='SGDID:S000000056',
                unique_id='SGD_0001',
                annotation_source_database='SGD',
                sgd_id='S000000056',
                ensembl_id='',
                guide_targets='sgd/1/23/TGCA/0/0',
            ),
        ])
        self.write_gpml([
            ('Ensembl', 'YAL060W', 'BDH1 Ensembl'),
            ('SGD', 'S000000056', 'BDH1 SGD'),
        ])
        analysis = PathwayAnalysis(self.config)

        pathway = analysis.process_all_pathways()[0]

        self.assertEqual(pathway['total_genes'], 1)
        self.assertEqual(pathway['matched_gene_count'], 0)
        self.assertEqual(pathway['unmatched_gene_count'], 1)
        self.assertEqual(pathway['genes_with_guides'], 0)
        self.assertEqual(pathway['pw_kfi_score'], 0.0)
        gene = pathway['unmatched_genes'][0]
        self.assertEqual(
            gene['reason'], 'ambiguous_canonical_annotation_identity'
        )
        self.assertEqual(gene['annotation_record_count'], 2)
        self.assertEqual(len(gene['annotation_identity_candidates']), 2)
        self.assertEqual(
            {xref['match_status'] for xref in gene['pathway_xrefs']},
            {'unmatched'},
        )

    def test_mixed_database_pathway_without_exact_crosslink_is_skipped(self):
        self.write_annotation([
            self.annotation_row(
                subject_id='CUSTOM:BDH1',
                annotation_source_database='Custom',
                sgd_id='',
            )
        ])
        self.write_gpml([
            ('Ensembl', 'YAL060W', 'BDH1 Ensembl'),
            ('SGD', 'S000000056', 'BDH1 SGD'),
        ])
        self.config['input']['sgd_db'] = str(self.temp_dir / 'missing.fasta')
        analysis = PathwayAnalysis(self.config)

        pathways = analysis.process_all_pathways()

        self.assertEqual(pathways, [])
        self.assertEqual(analysis.coverage_stats['skipped_pathways'], 1)
        self.assertEqual(len(analysis.pathway_gene_audit), 2)
        self.assertEqual(
            {
                row['pathway_canonicalization_status']
                for row in analysis.pathway_gene_audit
            },
            {'incomplete_cross_database_identity'},
        )
        ensembl_audit = next(
            row for row in analysis.pathway_gene_audit
            if row['pathway_gene_database'] == 'Ensembl'
        )
        self.assertFalse(ensembl_audit['canonicalization_complete'])

    def test_explicit_annotation_sgd_id_crosslinks_without_fasta_resolver(self):
        self.write_annotation([
            self.annotation_row(
                subject_id='CUSTOM:BDH1',
                annotation_source_database='Custom',
                sgd_id='S000000056',
                ensembl_id='YAL060W',
            )
        ])
        self.write_gpml([
            ('Ensembl', 'YAL060W', 'BDH1 Ensembl'),
            ('SGD', 'S000000056', 'BDH1 SGD'),
        ])
        self.config['input']['sgd_db'] = str(self.temp_dir / 'missing.fasta')
        analysis = PathwayAnalysis(self.config)

        pathway = analysis.process_all_pathways()[0]

        self.assertTrue(pathway['canonicalization_complete'])
        self.assertEqual(pathway['total_genes'], 1)
        self.assertEqual(pathway['matched_gene_count'], 1)
        self.assertEqual(pathway['pw_kfi_score'], 1.0)
        gene = pathway['matched_genes'][0]
        self.assertEqual(gene['canonical_subject_id'], 'SGDID:S000000056')
        self.assertIn(
            'exact_annotation_sgd_id', gene['canonicalization_methods']
        )

    def test_duplicate_unmatched_nodes_count_once_in_half_targetable_pathway(self):
        self.write_annotation([self.annotation_row()])
        self.write_gpml([
            ('SGD', 'S000000056', 'BDH1 guided first'),
            ('SGD', 'S000000056', 'BDH1 guided duplicate'),
            ('SGD', 'S000000057', 'unmatched first'),
            ('SGD', 'S000000057', 'unmatched duplicate'),
        ])
        analysis = PathwayAnalysis(self.config)

        pathway = analysis.process_all_pathways()[0]

        self.assertEqual(pathway['evaluable_gene_node_count'], 4)
        self.assertEqual(pathway['unique_xref_count'], 2)
        self.assertEqual(pathway['total_genes'], 2)
        self.assertEqual(pathway['matched_gene_count'], 1)
        self.assertEqual(pathway['unmatched_gene_count'], 1)
        self.assertEqual(pathway['genes_with_guides'], 1)
        self.assertEqual(pathway['pw_kfi_score'], 0.5)
        self.assertEqual(pathway['targetable_type'], 'Partially targetable')

    def test_duplicate_annotation_subject_aggregation_is_input_order_independent(self):
        rows = [
            self.annotation_row(
                query_id='query-z', evalue='1e-20', bit_score='100',
                guide_targets='g2/2/24/TGCA/0/0',
            ),
            self.annotation_row(
                query_id='query-a', evalue='1e-50', bit_score='200',
                guide_targets='g1/1/23/ACGT/0/0',
            ),
        ]
        self.write_gpml([('SGD', 'S000000056', 'BDH1')])

        self.write_annotation(rows)
        forward = PathwayAnalysis(self.config)
        forward_annotations = forward.load_annotation_results()
        forward_pathways = forward.process_all_pathways()

        self.write_annotation(list(reversed(rows)))
        reverse = PathwayAnalysis(self.config)
        reverse_annotations = reverse.load_annotation_results()
        reverse_pathways = reverse.process_all_pathways()

        forward_record = forward_annotations['by_subject_id'][
            'SGDID:S000000056'
        ]
        reverse_record = reverse_annotations['by_subject_id'][
            'SGDID:S000000056'
        ]
        self.assertEqual(forward_record, reverse_record)
        self.assertEqual(
            forward_record['annotation_query_ids'], ['query-a', 'query-z']
        )
        self.assertEqual(
            forward_record['guide_targets'],
            'g1/1/23/ACGT/0/0,g2/2/24/TGCA/0/0',
        )
        self.assertEqual(forward_record['query_id'], 'query-a')
        self.assertEqual(forward_pathways, reverse_pathways)
        self.assertEqual(
            forward.annotation_identifier_issues,
            reverse.annotation_identifier_issues,
        )

    def test_conflicting_subject_aliases_are_rejected_independently_of_row_order(self):
        rows = [
            self.annotation_row(query_id='query-a', ensembl_id='YAL060W'),
            self.annotation_row(query_id='query-b', ensembl_id='YBR160W'),
        ]
        observed = []
        for ordered_rows in (rows, list(reversed(rows))):
            self.write_annotation(ordered_rows)
            analysis = PathwayAnalysis(self.config)
            annotation_data = analysis.load_annotation_results()
            observed.append((annotation_data, analysis.annotation_identifier_issues))

        for annotation_data, issues in observed:
            self.assertIn(
                'SGDID:S000000056', annotation_data['ambiguous_subject_ids']
            )
            self.assertNotIn(
                'SGDID:S000000056', annotation_data['by_subject_id']
            )
            self.assertEqual(
                annotation_data['ambiguous_by_database']['Ensembl'],
                {'YAL060W', 'YBR160W'},
            )
            self.assertEqual(
                [issue['reason'] for issue in issues],
                ['ambiguous_annotation_subject_id'],
            )
        self.assertEqual(observed[0][1], observed[1][1])

    def test_run_saves_unmatched_audit_even_when_no_pathway_is_processed(self):
        self.write_annotation([
            self.annotation_row(
                subject_id='SGDID:S000000057',
                sgd_id='S000000057',
                ensembl_id='',
            )
        ])
        self.write_gpml([('Ensembl', 'YAL060W', 'BDH1')])
        self.config['pathway_analysis']['skip_low_coverage_pathways'] = True
        analysis = PathwayAnalysis(self.config)

        results = analysis.run()

        self.assertIsNotNone(results)
        self.assertEqual(results['processed_pathways'], [])
        audit_path = Path(results['output_files']['gene_match_audit_path'])
        self.assertTrue(audit_path.exists())
        with audit_path.open(encoding='utf-8', newline='') as handle:
            audit_rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(audit_rows), 1)
        self.assertEqual(audit_rows[0]['pathway_gene_database'], 'Ensembl')
        self.assertEqual(audit_rows[0]['status'], 'unmatched')
        self.assertEqual(audit_rows[0]['reason'], 'annotation_not_found')
        for output_path in results['output_files'].values():
            self.assertTrue(Path(output_path).is_file())
            self.assertGreater(Path(output_path).stat().st_size, 0)

        report_results = PathwayReportGenerator(self.config).run()

        self.assertNotIn('database_path', report_results)
        self.assertTrue(
            Path(report_results['report_sources_manifest']).is_file()
        )
        self.assertTrue(Path(report_results['summary_json_path']).is_file())
        self.assertTrue(Path(report_results['html_report_path']).is_file())
        self.assertEqual(report_results['statistics']['total_pathways'], 0)

    @unittest.skipUnless(
        BUNDLED_GPML_DIR.is_dir(),
        'Bundled WikiPathways GPML data is not available.',
    )
    def test_bundled_gpml_retains_all_strict_ensembl_nodes_and_only_pathways(self):
        self.config['pathway_analysis']['wikipathways_gpml_dir'] = str(
            BUNDLED_GPML_DIR
        )
        analysis = PathwayAnalysis(self.config)
        pathways = [
            analysis.parse_gpml_file(path)
            for path in analysis.get_gpml_files()
        ]
        pathways = [pathway for pathway in pathways if pathway is not None]

        ensembl_nodes = sum(
            gene['database'] == 'Ensembl'
            for pathway in pathways
            for gene in pathway['genes']
        )
        ensembl_only_pathways = {
            pathway['pathway_id']
            for pathway in pathways
            if any(gene['database'] == 'Ensembl' for gene in pathway['genes'])
            and not any(gene['database'] == 'SGD' for gene in pathway['genes'])
        }

        self.assertEqual(ensembl_nodes, 254)
        self.assertEqual(
            ensembl_only_pathways,
            {'WP5201', 'WP14', 'WP490', 'WP112', 'WP515', 'WP102', 'WP369'},
        )
        wp14 = next(
            pathway for pathway in pathways if pathway['pathway_id'] == 'WP14'
        )
        self.assertEqual(wp14['total_genes'], 1)
        self.assertEqual(
            {row['reason'] for row in wp14['unmatchable_genes']},
            {'invalid_identifier_format'},
        )

    @unittest.skipUnless(
        BUNDLED_GPML_DIR.is_dir() and BUNDLED_SGD_PATH.is_file(),
        'Bundled WikiPathways and SGD data are not available.',
    )
    def test_bundled_gpml_exact_gene_count_contract(self):
        self.config['input']['sgd_db'] = str(BUNDLED_SGD_PATH)
        self.config['pathway_analysis']['wikipathways_gpml_dir'] = str(
            BUNDLED_GPML_DIR
        )
        analysis = PathwayAnalysis(self.config)
        resolver = analysis._load_sgd_identifier_crosswalk()
        annotation_data = {'cross_references': {'Ensembl': resolver}}

        raw_nodes = 0
        evaluable_nodes = 0
        unique_xrefs = 0
        canonical_genes = 0
        for gpml_path in analysis.get_gpml_files():
            pathway = analysis.parse_gpml_file(gpml_path)
            self.assertIsNotNone(pathway)
            raw_nodes += pathway['total_gene_nodes']
            evaluable_nodes += pathway['evaluable_gene_node_count']
            unique_xrefs += pathway['unique_xref_count']
            keys = set()
            for gene in pathway['genes']:
                subject_id, _, complete = (
                    analysis._canonical_subject_for_pathway_gene(
                        gene, annotation_data
                    )
                )
                self.assertTrue(complete)
                keys.add(analysis._subject_gene_key(subject_id))
            canonical_genes += len(keys)

        self.assertEqual(raw_nodes, 1602)
        self.assertEqual(evaluable_nodes, 1242)
        self.assertEqual(unique_xrefs, 1132)
        self.assertEqual(canonical_genes, 1131)


if __name__ == '__main__':
    unittest.main()
