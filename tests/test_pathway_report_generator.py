import csv
import json
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from src.modules.pathway_report_generator import PathwayReportGenerator


class PathwayReportIdentifierTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir_context = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temp_dir_context.name)
        self.intermediate_dir = self.temp_dir / 'intermediate'
        self.output_dir = self.temp_dir / 'output'
        self.intermediate_dir.mkdir()
        self.generator = PathwayReportGenerator({
            'input': {
                'genome_assembly': 'fixture_genome.fa',
                'transcriptome': str(
                    self.intermediate_dir / 'transcriptome.fa'
                ),
                'protein_db_for_goanalysis': 'fixture_protein_db.fa',
            },
            'crispr': {
                'pam_pattern': 'NAG',
                'max_mismatch': 3,
                'threads': 7,
            },
            'busco': {'lineage': 'fixture_lineage'},
            'database': {'force_rebuild': False},
            'report': {
                'generate_summary': False,
                'generate_html': True,
            },
            'output': {
                'intermediate': str(self.intermediate_dir),
                'base_dir': str(self.output_dir),
            },
            'pathway_report': {},
        })

    def tearDown(self):
        self.temp_dir_context.cleanup()

    @staticmethod
    def matched_gene(**overrides):
        gene = {
            'pathway_gene_database': 'Ensembl',
            'pathway_gene_id': 'YAL060W',
            'pathway_gene_name': 'BDH1',
            'gene_key': 'Ensembl::YAL060W',
            'pw_kfi_gene_key': 'SGD::SGDID:S000000056',
            'pathway_xrefs': repr([{
                'database': 'Ensembl',
                'gene_id': 'YAL060W',
                'gene_key': 'Ensembl::YAL060W',
                'annotation_key': 'ensembl_id',
                'pathway_gene_name': 'BDH1',
                'pathway_gene_names': ['BDH1'],
                'graph_ids': ['n1'],
                'node_occurrence_count': 1,
                'match_status': 'matched',
                'reason': '',
                'match_method': 'exact_sgd_fasta_cross_reference',
                'mapping_source': '/data/sgd.fasta',
                'matched_annotation_key': 'subject_id',
                'matched_annotation_id': 'SGDID:S000000056',
            }]),
            'source_databases': 'Ensembl',
            'node_occurrence_count': 1,
            'unique_xref_count': 1,
            'canonical_subject_id': 'SGDID:S000000056',
            'canonicalization_methods': 'sgd_fasta_cross_reference',
            'canonicalization_complete': True,
            'annotation_source_database': 'SGD',
            'matched_annotation_key': 'subject_id',
            'matched_annotation_id': 'SGDID:S000000056',
            'resolved_annotation_subject_id': 'SGDID:S000000056',
            'match_method': 'exact_sgd_fasta_cross_reference',
            'mapping_source': '/data/sgd.fasta',
            'sgd_id': 'S000000056',
            'ensembl_id': 'YAL060W',
            'unique_id': 'Single_0000000001',
            'gene_name': 'BDH1',
            'description': 'NAD-dependent protein',
            'guide_targets': 'guide-1/1/23/ACGT/0/0/100/+',
            'has_guides': True,
            'guide_count': 1,
        }
        gene.update(overrides)
        return gene

    def write_canonical_artifacts(
        self, *, gene_overrides=None, pathway_overrides=None,
        mapping_overrides=None,
    ):
        gene = self.matched_gene(**(gene_overrides or {}))
        pathway = {
            'pathway_id': 'WP9999',
            'pathway_name': 'Strict pathway',
            'pathway_organism': 'Saccharomyces cerevisiae',
            'total_gene_nodes': 1,
            'evaluable_gene_node_count': 1,
            'unique_xref_count': 1,
            'total_genes': 1,
            'gene_count_basis': 'exact_canonical_gene_units',
            'canonicalization_complete': True,
            'canonicalization_status': 'complete',
            'unmatchable_gene_count': 0,
            'unmatchable_genes': repr([]),
            'matched_genes': repr([gene]),
            'unmatched_genes': repr([]),
            'matched_gene_count': 1,
            'unmatched_gene_count': 0,
            'genes_with_guides': 1,
            'pw_kfi_score': 1.0,
            'targetable_type': 'Fully targetable',
            'coverage_rate': 1.0,
            'match_counts_by_database': repr({'Ensembl': 1}),
            'unmatched_counts_by_database': repr({}),
            'gpml_file_path': '/data/WP9999.gpml',
        }
        pathway.update(pathway_overrides or {})
        pathway_path = self.intermediate_dir / 'pathway_kfi_table.csv'
        with pathway_path.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(pathway))
            writer.writeheader()
            writer.writerow(pathway)

        def serialized_list(value):
            if isinstance(value, str):
                return value
            return ','.join(value)

        mapping = {
            'pathway_id': pathway['pathway_id'],
            'pathway_name': pathway['pathway_name'],
            'pathway_gene_database': gene['pathway_gene_database'],
            'pathway_gene_id': gene['pathway_gene_id'],
            'pathway_gene_name': gene['pathway_gene_name'],
            'gene_key': gene['gene_key'],
            'pw_kfi_gene_key': gene['pw_kfi_gene_key'],
            'pathway_xrefs': repr(
                self.generator._pathway_xrefs(gene)
            ),
            'source_databases': serialized_list(gene['source_databases']),
            'node_occurrence_count': gene['node_occurrence_count'],
            'unique_xref_count': gene['unique_xref_count'],
            'canonical_subject_id': gene['canonical_subject_id'],
            'canonicalization_methods': serialized_list(
                gene['canonicalization_methods']
            ),
            'canonicalization_complete': gene['canonicalization_complete'],
            'annotation_source_database': gene[
                'annotation_source_database'
            ],
            'matched_annotation_key': gene['matched_annotation_key'],
            'matched_annotation_id': gene['matched_annotation_id'],
            'resolved_annotation_subject_id': gene[
                'resolved_annotation_subject_id'
            ],
            'match_method': gene['match_method'],
            'mapping_source': gene['mapping_source'],
            'sgd_id': gene['sgd_id'],
            'ensembl_id': gene['ensembl_id'],
            'gene_sgd_id': gene['sgd_id'],
            'unique_id': gene['unique_id'],
            'gene_name': gene['gene_name'],
            'has_guides': gene['has_guides'],
            'guide_count': gene['guide_count'],
        }
        mapping.update(mapping_overrides or {})
        mapping_path = self.intermediate_dir / 'pathway_gene_mapping.tsv'
        with mapping_path.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(mapping), delimiter='\t'
            )
            writer.writeheader()
            writer.writerow(mapping)
        return pathway, mapping

    @staticmethod
    def write_gpml(path, pathway_name='Strict pathway'):
        path.write_text(
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Pathway xmlns="http://pathvisio.org/GPML/2013a" '
            f'Name="{pathway_name}" Organism="Saccharomyces cerevisiae">'
            '<DataNode TextLabel="BDH1" GraphId="n1" Type="GeneProduct">'
            '<Graphics CenterX="100" CenterY="100" Width="80" Height="20" />'
            '<Xref Database="Ensembl" ID="YAL060W" />'
            '</DataNode></Pathway>',
            encoding='utf-8',
        )
        return path

    @staticmethod
    def embedded_pathway_data(report_path):
        content = Path(report_path).read_text(encoding='utf-8')
        marker = '<script type="application/json" id="pathway-data">'
        payload = content.split(marker, 1)[1].split('</script>', 1)[0]
        return json.loads(payload)

    @staticmethod
    def embedded_report_data(report_path):
        content = Path(report_path).read_text(encoding='utf-8')
        marker = '<script type="application/json" id="twheel-data">'
        payload = content.split(marker, 1)[1].split('</script>', 1)[0]
        return json.loads(payload)

    def write_common_report_artifacts(self):
        (self.intermediate_dir / 'transcriptome.fa').write_text(
            '>fixture_transcript\nACGTACGT\n', encoding='utf-8'
        )
        (self.intermediate_dir / 'query.pep').write_text(
            '>fixture_protein\nMPEPTIDE\n', encoding='utf-8'
        )

        annotation_header = '\t'.join([
            'query_id', 'subject_id', 'identity', 'alignment_length',
            'evalue', 'bit_score', 'description', 'gene_name',
            'classification', 'unique_id',
        ])
        annotation_row = '\t'.join([
            'fixture_query', 'fixture_subject', '99.0', '100', '1e-20',
            '200', 'Fixture annotation description', 'FIXTURE_GENE',
            'Single', 'Single_fixture',
        ])
        (self.intermediate_dir / 'annotation_classified.txt').write_text(
            annotation_header + '\n' + annotation_row + '\n',
            encoding='utf-8',
        )
        (
            self.intermediate_dir / 'single_genes_targets_specific_clean.bed'
        ).write_text(
            'Single_fixture\t10\t33\n', encoding='utf-8'
        )
        (
            self.intermediate_dir / 'multi_genes_targets_specific_clean.bed'
        ).write_text('', encoding='utf-8')

        final_header = annotation_header + '\tguide_targets'
        final_row = annotation_row + (
            '\tSingle_fixture/10/33/FIXTURE_GUIDE_NAG/1/0/0.987/+'
        )
        (self.intermediate_dir / 'annotation_with_guides_final.tsv').write_text(
            final_header + '\n' + final_row + '\n', encoding='utf-8'
        )

    def test_svg_uses_exact_database_and_identifier_composite_key(self):
        pathway_data = {
            'board_width': 400,
            'board_height': 300,
            'edges': [],
            'nodes': [
                {
                    'type': 'GeneProduct',
                    'label': 'BDH1',
                    'xref_id': 'YAL060W',
                    'database': 'Ensembl',
                    'x': 100,
                    'y': 100,
                    'width': 80,
                    'height': 20,
                },
                {
                    'type': 'GeneProduct',
                    'label': 'wrong database',
                    'xref_id': 'YAL060W',
                    'database': 'SGD',
                    'x': 200,
                    'y': 100,
                    'width': 80,
                    'height': 20,
                },
            ],
        }
        status_map = {
            'Ensembl::YAL060W': {
                'has_guides': True,
                'sgd_id': 'S000000056',
            }
        }

        svg = self.generator.generate_pathway_svg(pathway_data, status_map)

        self.assertIn('gene-node gene-with-guides', svg)
        self.assertIn('data-gene-key="Ensembl::YAL060W"', svg)
        self.assertIn('data-gene-database="Ensembl"', svg)
        self.assertIn('gene-node gene-unknown', svg)
        self.assertNotIn('guide.gene_name', svg)

    def test_guide_csv_data_preserves_pathway_and_annotation_provenance(self):
        pathway = {
            'pathway_id': 'WP9999',
            'pathway_name': 'Strict pathway',
            'matched_genes': repr([self.matched_gene()]),
        }

        result = self.generator.generate_pathway_guide_csv_data([pathway])

        guide = result['WP9999']['guides'][0]
        self.assertEqual(guide['gene_key'], 'Ensembl::YAL060W')
        self.assertEqual(guide['pathway_gene_database'], 'Ensembl')
        self.assertEqual(guide['pathway_gene_id'], 'YAL060W')
        self.assertEqual(guide['annotation_source_database'], 'SGD')
        self.assertEqual(guide['matched_annotation_id'], 'SGDID:S000000056')
        self.assertEqual(guide['sgd_id'], 'S000000056')
        self.assertEqual(guide['ensembl_id'], 'YAL060W')
        self.assertEqual(guide['MM0'], '0')
        self.assertEqual(guide['MM1-2'], '0')
        self.assertNotIn('mm0', guide)
        self.assertNotIn('mm1_2', guide)

    def test_diagram_exposes_exact_unmatched_reason_without_coloring_as_match(self):
        gpml_path = self.temp_dir / 'WP9998.gpml'
        gpml_path.write_text(
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Pathway xmlns="http://pathvisio.org/GPML/2013a" '
            'Name="Unmatched pathway" Organism="Saccharomyces cerevisiae">'
            '<DataNode TextLabel="BDH1" GraphId="n1" Type="GeneProduct">'
            '<Graphics CenterX="100" CenterY="100" Width="80" Height="20" />'
            '<Xref Database=" Ensembl " ID=" YAL060W " />'
            '</DataNode></Pathway>',
            encoding='utf-8',
        )
        pathways = [{
            'pathway_id': 'WP9998',
            'pathway_name': 'Unmatched pathway',
            'gpml_file_path': str(gpml_path),
            'matched_genes': '[]',
            'unmatched_genes': repr([{
                'database': 'Ensembl',
                'gene_id': 'YAL060W',
                'gene_name': 'BDH1',
                'reason': 'annotation_not_found',
            }]),
            'unmatchable_genes': '[]',
        }]

        diagram = self.generator.get_pathway_diagram_data('WP9998', pathways)

        status = diagram['guide_status_map']['Ensembl::YAL060W']
        self.assertEqual(status['match_status'], 'unmatched')
        self.assertEqual(status['failure_reason'], 'annotation_not_found')
        self.assertIn(
            'class="gene-node gene-unknown" data-gene-database="Ensembl"',
            diagram['svg_content'],
        )

    def test_diagram_colors_duplicate_nodes_but_counts_one_pw_kfi_gene(self):
        gpml_path = self.temp_dir / 'WP9997.gpml'
        gpml_path.write_text(
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Pathway xmlns="http://pathvisio.org/GPML/2013a" '
            'Name="Duplicate nodes" Organism="Saccharomyces cerevisiae">'
            '<DataNode TextLabel="BDH1-a" GraphId="n1" Type="GeneProduct">'
            '<Graphics CenterX="100" CenterY="100" Width="80" Height="20" />'
            '<Xref Database="Ensembl" ID="YAL060W" /></DataNode>'
            '<DataNode TextLabel="BDH1-b" GraphId="n2" Type="GeneProduct">'
            '<Graphics CenterX="200" CenterY="100" Width="80" Height="20" />'
            '<Xref Database="Ensembl" ID="YAL060W" />'
            '</DataNode>'
            '<DataNode TextLabel="unsupported" GraphId="n3" Type="GeneProduct">'
            '<Graphics CenterX="300" CenterY="100" Width="80" Height="20" />'
            '<Xref Database="Entrez Gene" ID="123" />'
            '</DataNode></Pathway>',
            encoding='utf-8',
        )
        pathways = [{
            'pathway_id': 'WP9997',
            'pathway_name': 'Duplicate nodes',
            'gpml_file_path': str(gpml_path),
            'matched_genes': repr([self.matched_gene(
                pw_kfi_gene_key='SGD::SGDID:S000000056',
                pathway_xrefs=[{
                    'database': 'Ensembl',
                    'gene_id': 'YAL060W',
                    'gene_key': 'Ensembl::YAL060W',
                    'node_occurrence_count': 2,
                    'match_status': 'matched',
                }],
                node_occurrence_count=2,
                unique_xref_count=1,
            )]),
            'unmatched_genes': '[]',
            'unmatchable_genes': '[]',
            'total_genes': '1',
            'evaluable_gene_node_count': '2',
            'unique_xref_count': '1',
            'genes_with_guides': '1',
        }]

        diagram = self.generator.get_pathway_diagram_data('WP9997', pathways)

        self.assertEqual(diagram['total_genes'], 1)
        self.assertEqual(diagram['evaluable_gene_nodes'], 2)
        self.assertEqual(diagram['total_gene_nodes'], 3)
        self.assertEqual(diagram['unique_xref_count'], 1)
        self.assertEqual(diagram['genes_with_guides'], 1)
        self.assertEqual(diagram['gene_nodes_with_guides'], 2)
        self.assertEqual(
            diagram['guide_status_map']['Ensembl::YAL060W'][
                'node_occurrence_count'
            ],
            2,
        )

    def test_cross_database_alias_colors_both_nodes_and_exports_guide_once(self):
        gpml_path = self.temp_dir / 'WP9994.gpml'
        gpml_path.write_text(
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Pathway xmlns="http://pathvisio.org/GPML/2013a" '
            'Name="Cross database alias" Organism="Saccharomyces cerevisiae">'
            '<DataNode TextLabel="BDH1 Ensembl" GraphId="n1" Type="GeneProduct">'
            '<Graphics CenterX="100" CenterY="100" Width="80" Height="20" />'
            '<Xref Database="Ensembl" ID="YAL060W" /></DataNode>'
            '<DataNode TextLabel="BDH1 SGD" GraphId="n2" Type="GeneProduct">'
            '<Graphics CenterX="200" CenterY="100" Width="80" Height="20" />'
            '<Xref Database="SGD" ID="S000000056" /></DataNode>'
            '</Pathway>',
            encoding='utf-8',
        )
        matched_gene = self.matched_gene(
            pw_kfi_gene_key='SGD::SGDID:S000000056',
            pathway_xrefs=[
                {
                    'database': 'Ensembl',
                    'gene_id': 'YAL060W',
                    'gene_key': 'Ensembl::YAL060W',
                    'match_status': 'matched',
                },
                {
                    'database': 'SGD',
                    'gene_id': 'S000000056',
                    'gene_key': 'SGD::S000000056',
                    'match_status': 'matched',
                },
            ],
        )
        pathways = [{
            'pathway_id': 'WP9994',
            'pathway_name': 'Cross database alias',
            'gpml_file_path': str(gpml_path),
            'matched_genes': repr([matched_gene]),
            'unmatched_genes': '[]',
            'unmatchable_genes': '[]',
            'total_genes': '1',
            'evaluable_gene_node_count': '2',
            'unique_xref_count': '2',
            'genes_with_guides': '1',
        }]

        diagram = self.generator.get_pathway_diagram_data('WP9994', pathways)
        guides = self.generator.generate_pathway_guide_csv_data(pathways)

        self.assertEqual(diagram['total_genes'], 1)
        self.assertEqual(diagram['genes_with_guides'], 1)
        self.assertEqual(diagram['gene_nodes_with_guides'], 2)
        self.assertEqual(
            diagram['guide_status_map']['Ensembl::YAL060W'][
                'pw_kfi_gene_key'
            ],
            'SGD::SGDID:S000000056',
        )
        self.assertEqual(
            diagram['guide_status_map']['SGD::S000000056'][
                'pw_kfi_gene_key'
            ],
            'SGD::SGDID:S000000056',
        )
        self.assertEqual(len(guides['WP9994']['guides']), 1)
        self.assertEqual(
            guides['WP9994']['guides'][0]['source_xref_keys'],
            ['Ensembl::YAL060W', 'SGD::S000000056'],
        )

    def test_legacy_node_count_pathway_table_is_rejected(self):
        kfi_path = self.intermediate_dir / 'pathway_kfi_table.csv'
        kfi_path.write_text(
            'pathway_id,pathway_name,total_genes,genes_with_guides,pw_kfi_score\n'
            'WP1,Legacy,2,2,1.0\n',
            encoding='utf-8',
        )
        (self.intermediate_dir / 'pathway_gene_mapping.tsv').write_text(
            '', encoding='utf-8'
        )

        with self.assertRaisesRegex(
            ValueError, 'Legacy or incomplete PW-KFI table'
        ):
            self.generator.load_pathway_analysis_results()

    def test_canonical_artifacts_reconcile_before_report_loading(self):
        self.write_canonical_artifacts()

        pathways, pathway_genes = (
            self.generator.load_pathway_analysis_results()
        )

        self.assertEqual(len(pathways), 1)
        self.assertEqual(
            pathway_genes['WP9999'][0]['pw_kfi_gene_key'],
            'SGD::SGDID:S000000056',
        )

    def test_single_database_unit_does_not_require_cross_database_resolution(self):
        self.write_canonical_artifacts(
            gene_overrides={'canonicalization_complete': False}
        )

        pathways, pathway_genes = (
            self.generator.load_pathway_analysis_results()
        )

        self.assertEqual(len(pathways), 1)
        self.assertEqual(len(pathway_genes['WP9999']), 1)

    def test_cross_database_units_require_exact_canonical_resolution(self):
        ensembl_xref = self.generator._pathway_xrefs(
            self.matched_gene()
        )[0]
        self.write_canonical_artifacts(
            gene_overrides={
                'pathway_xrefs': repr([
                    ensembl_xref,
                    {
                        'database': 'SGD',
                        'gene_id': 'S000000056',
                        'gene_key': 'SGD::S000000056',
                        'node_occurrence_count': 1,
                    },
                ]),
                'source_databases': 'Ensembl,SGD',
                'node_occurrence_count': 2,
                'unique_xref_count': 2,
                'canonicalization_complete': False,
            },
            pathway_overrides={
                'total_gene_nodes': 2,
                'evaluable_gene_node_count': 2,
                'unique_xref_count': 2,
            },
        )

        with self.assertRaisesRegex(
            ValueError, 'Cross-database canonicalization does not reconcile'
        ):
            self.generator.load_pathway_analysis_results()

    def test_mapping_payload_mismatch_is_rejected_with_same_canonical_key(self):
        self.write_canonical_artifacts(mapping_overrides={
            'pathway_xrefs': repr([{
                'database': 'SGD',
                'gene_id': 'S999999999',
                'gene_key': 'SGD::S999999999',
                'node_occurrence_count': 999,
            }]),
            'node_occurrence_count': 999,
            'canonical_subject_id': 'SGDID:S999999999',
            'canonicalization_complete': False,
        })

        with self.assertRaisesRegex(
            ValueError, 'Pathway mapping payload does not reconcile'
        ):
            self.generator.load_pathway_analysis_results()

    def test_guide_count_mismatch_in_canonical_gene_is_rejected(self):
        self.write_canonical_artifacts(
            gene_overrides={'guide_count': 0}
        )

        with self.assertRaisesRegex(
            ValueError, 'Guide counts do not reconcile'
        ):
            self.generator.load_pathway_analysis_results()

    def test_malformed_serialized_gene_data_fails_instead_of_rendering_unknown(self):
        gpml_path = self.temp_dir / 'WP9996.gpml'
        gpml_path.write_text(
            '<?xml version="1.0" encoding="UTF-8"?>'
            '<Pathway xmlns="http://pathvisio.org/GPML/2013a" Name="Broken">'
            '<DataNode TextLabel="BDH1" GraphId="n1" Type="GeneProduct">'
            '<Xref Database="Ensembl" ID="YAL060W" />'
            '</DataNode></Pathway>',
            encoding='utf-8',
        )
        pathways = [{
            'pathway_id': 'WP9996',
            'pathway_name': 'Broken',
            'gpml_file_path': str(gpml_path),
            'matched_genes': 'not valid serialized data',
            'unmatched_genes': '[]',
            'unmatchable_genes': '[]',
            'targetable_type': 'Not targetable',
        }]

        with self.assertRaisesRegex(ValueError, 'Cannot deserialize'):
            self.generator.get_pathway_diagram_data('WP9996', pathways)
        with self.assertRaisesRegex(ValueError, 'Cannot deserialize matched genes'):
            self.generator.generate_pathway_guide_csv_data(pathways)

        summary = self.generator.generate_pathway_json_summary(pathways)
        self.assertEqual(summary['statistics']['svg_generation']['failed'], 1)
        self.assertIn(
            'Cannot deserialize',
            summary['statistics']['svg_generation']['failure_reasons']['WP9996'],
        )

    def test_run_uses_staged_sources_after_intermediate_and_gpml_are_removed(self):
        gpml_path = self.write_gpml(self.temp_dir / 'WP9999.gpml')
        self.write_canonical_artifacts(
            pathway_overrides={'gpml_file_path': str(gpml_path)}
        )
        database_path = Path(self.generator.web_data_dir) / 'twheel_data.db'
        database_path.write_bytes(b'pre-existing database sentinel')

        first_results = self.generator.run()
        first_embedded = self.embedded_pathway_data(
            first_results['html_report_path']
        )

        self.assertNotIn('database_path', first_results)
        self.assertEqual(
            database_path.read_bytes(), b'pre-existing database sentinel'
        )
        manifest_path = Path(first_results['report_sources_manifest'])
        snapshot = self.generator.report_source_store.load_namespace('pathway')
        self.assertIn('pathway_kfi_table', snapshot)
        self.assertIn('pathway_gene_mapping', snapshot)
        self.assertIn('gpml.WP9999', snapshot)

        shutil.rmtree(self.intermediate_dir)
        gpml_path.unlink()
        second_results = PathwayReportGenerator(self.generator.config).run()
        second_embedded = self.embedded_pathway_data(
            second_results['html_report_path']
        )
        summary = json.loads(
            Path(second_results['summary_json_path']).read_text(encoding='utf-8')
        )

        first_embedded['summary'].pop('timestamp', None)
        second_embedded['summary'].pop('timestamp', None)
        self.assertEqual(second_embedded, first_embedded)
        self.assertEqual(summary['statistics']['svg_generation']['success'], 1)
        self.assertTrue(Path(second_results['html_report_path']).is_file())
        self.assertTrue(manifest_path.is_file())
        self.assertFalse(self.intermediate_dir.exists())
        self.assertEqual(
            database_path.read_bytes(), b'pre-existing database sentinel'
        )

    def test_invalid_current_pair_preserves_previous_snapshot(self):
        gpml_path = self.write_gpml(self.temp_dir / 'WP9999.gpml')
        self.write_canonical_artifacts(
            pathway_overrides={'gpml_file_path': str(gpml_path)}
        )
        self.generator.run()
        manifest_path = self.generator.report_source_store.manifest_path
        manifest_before = manifest_path.read_bytes()
        snapshot_before = {
            key: Path(path).read_bytes()
            for key, path in self.generator.report_source_store.load_namespace(
                'pathway'
            ).items()
        }
        (self.intermediate_dir / 'pathway_kfi_table.csv').write_text(
            'pathway_id,pathway_name\nWP9999,broken\n',
            encoding='utf-8',
        )

        with self.assertRaisesRegex(
            ValueError, 'Legacy or incomplete PW-KFI table'
        ):
            self.generator.load_pathway_analysis_results()

        self.assertEqual(manifest_path.read_bytes(), manifest_before)
        snapshot_after = self.generator.report_source_store.load_namespace(
            'pathway'
        )
        self.assertEqual(
            {key: Path(path).read_bytes() for key, path in snapshot_after.items()},
            snapshot_before,
        )
        shutil.rmtree(self.intermediate_dir)
        regenerated = PathwayReportGenerator(self.generator.config).run()
        summary = json.loads(
            Path(regenerated['summary_json_path']).read_text(encoding='utf-8')
        )
        self.assertEqual(summary['statistics']['svg_generation']['success'], 1)

    def test_new_tables_with_missing_gpml_do_not_reuse_old_snapshot(self):
        gpml_path = self.write_gpml(self.temp_dir / 'WP9999.gpml')
        self.write_canonical_artifacts(
            pathway_overrides={'gpml_file_path': str(gpml_path)}
        )
        self.generator.run()
        first_snapshot = self.generator.report_source_store.load_namespace(
            'pathway'
        )
        self.assertIn('gpml.WP9999', first_snapshot)

        self.write_canonical_artifacts(pathway_overrides={
            'gpml_file_path': str(self.temp_dir / 'missing-WP9999.gpml')
        })
        results = self.generator.run()
        summary = json.loads(
            Path(results['summary_json_path']).read_text(encoding='utf-8')
        )
        second_snapshot = self.generator.report_source_store.load_namespace(
            'pathway'
        )

        self.assertNotIn('gpml.WP9999', second_snapshot)
        self.assertEqual(summary['statistics']['svg_generation']['success'], 0)
        self.assertEqual(summary['statistics']['svg_generation']['no_gpml'], 1)

    def test_existing_but_empty_or_partial_intermediate_never_falls_back(self):
        self.write_canonical_artifacts()
        self.generator.load_pathway_analysis_results()
        mapping_path = self.intermediate_dir / 'pathway_gene_mapping.tsv'
        kfi_path = self.intermediate_dir / 'pathway_kfi_table.csv'

        mapping_path.unlink()
        with self.assertRaisesRegex(
            FileNotFoundError,
            'Incomplete pathway report source snapshot in intermediate',
        ):
            self.generator.load_pathway_analysis_results()

        kfi_path.unlink()
        with self.assertRaisesRegex(
            FileNotFoundError,
            'Pathway report sources are absent from intermediate',
        ):
            self.generator.load_pathway_analysis_results()

    def test_gpml_snapshot_names_resist_sanitization_collisions(self):
        first_key, first_path = self.generator._gpml_source_entry('A/B')
        second_key, second_path = self.generator._gpml_source_entry('A_B')

        self.assertNotEqual(first_key, second_key)
        self.assertNotEqual(first_path, second_path)
        self.assertTrue(first_path.startswith('pathway/gpml/A_B-'))
        self.assertTrue(second_path.startswith('pathway/gpml/A_B-'))

    def test_duplicate_canonical_mapping_rows_are_rejected(self):
        self.write_canonical_artifacts()
        mapping_path = self.intermediate_dir / 'pathway_gene_mapping.tsv'
        lines = mapping_path.read_text(encoding='utf-8').splitlines()
        mapping_path.write_text(
            '\n'.join([lines[0], lines[1], lines[1]]) + '\n',
            encoding='utf-8',
        )

        with self.assertRaisesRegex(
            ValueError, 'Pathway mapping does not reconcile'
        ):
            self.generator.load_pathway_analysis_results()

    def test_unparseable_matched_gene_count_is_rejected_by_loader(self):
        self.write_canonical_artifacts(pathway_overrides={
            'matched_gene_count': 'invalid'
        })

        with self.assertRaisesRegex(
            ValueError, 'Invalid matched_gene_count for WP9999'
        ):
            self.generator.load_pathway_analysis_results()

    def test_html_sidebar_and_sections_follow_the_go_report_contract(self):
        html_content = self.generator.create_pathway_html_template(
            [],
            {'pathway_guide_data': {}, 'pathway_svgs': {}, 'statistics': {}},
            {},
        )

        nav_items = re.findall(
            r'<li>\s*<a href="#([^"]+)" '
            r'class="nav-link(?: active)?">([^<]+)</a>\s*</li>',
            html_content,
        )
        expected_items = [
            ('general-stats', 'General Statistics'),
            ('gene-prediction', 'Gene Prediction'),
            ('gene-annotation', 'Gene Annotation'),
            ('guide-rna', 'Guide RNA Design'),
            ('pw-kfi-assessment', 'PW-KFI Assessment'),
            ('configuration', 'Analysis Configuration'),
        ]
        self.assertEqual(nav_items, expected_items)
        self.assertEqual(
            re.findall(
                r'<section id="([^"]+)" '
                r'class="content-section(?: [^"]*)?">',
                html_content,
            ),
            [section_id for section_id, _ in expected_items],
        )
        self.assertIn(
            '<a href="#general-stats" class="nav-link active">',
            html_content,
        )
        for _, label in expected_items:
            with self.subTest(section_heading=label):
                self.assertEqual(
                    html_content.count(f'<h2>{label}</h2>'), 1
                )

    def test_pw_kfi_content_is_grouped_in_one_assessment_section(self):
        html_content = self.generator.create_pathway_html_template(
            [],
            {'pathway_guide_data': {}, 'pathway_svgs': {}, 'statistics': {}},
            {},
        )

        assessment_start = html_content.index(
            '<section id="pw-kfi-assessment"'
        )
        assessment_end = html_content.index(
            '<section id="configuration" class="content-section">',
            assessment_start,
        )
        assessment_html = html_content[assessment_start:assessment_end]
        outside_assessment = (
            html_content[:assessment_start] + html_content[assessment_end:]
        )

        self.assertEqual(
            assessment_html.count(
                '<h3>Pathway Targetability Distribution</h3>'
            ),
            1,
        )
        self.assertEqual(
            assessment_html.count('<h3>Pathway Analysis Results</h3>'),
            1,
        )
        self.assertNotIn(
            '<h3>Pathway Targetability Distribution</h3>', outside_assessment
        )
        self.assertNotIn(
            '<h3>Pathway Analysis Results</h3>', outside_assessment
        )
        self.assertEqual(
            html_content.count('<h2>PW-KFI Assessment</h2>'), 1
        )
        self.assertNotIn('<section id="summary-stats"', html_content)
        self.assertNotIn('<section id="pathway-analysis"', html_content)
        self.assertNotIn('<h2>Summary Statistics</h2>', html_content)
        self.assertNotIn('<h2>Pathway Analysis</h2>', html_content)

    def test_common_sections_embed_run_data_and_analysis_configuration(self):
        self.write_canonical_artifacts()
        self.write_common_report_artifacts()

        results = self.generator.run()
        html_content = Path(results['html_report_path']).read_text(
            encoding='utf-8'
        )
        embedded_data = self.embedded_report_data(
            results['html_report_path']
        )

        stats = embedded_data['statistics']
        self.assertEqual(stats['gene_prediction']['transcript_count'], 1)
        self.assertEqual(stats['gene_prediction']['protein_count'], 1)
        self.assertEqual(
            stats['gene_prediction']['proteins'][0]['id'],
            'fixture_protein',
        )
        self.assertEqual(stats['annotation']['total_genes'], 1)
        self.assertEqual(stats['annotation']['single_genes'], 1)
        self.assertEqual(stats['guide_design']['total_guides'], 1)

        self.assertEqual(len(embedded_data['gene_annotation']), 1)
        self.assertEqual(
            embedded_data['gene_annotation'][0]['gene_name'],
            'FIXTURE_GENE',
        )
        self.assertEqual(
            embedded_data['gene_annotation'][0]['description'],
            'Fixture annotation description',
        )
        self.assertEqual(len(embedded_data['guide_rna']), 1)
        self.assertEqual(
            embedded_data['guide_rna'][0]['target_sequence'],
            'FIXTURE_GUIDE_NAG',
        )
        self.assertEqual(embedded_data['guide_rna'][0]['mm0'], 1)
        self.assertEqual(embedded_data['guide_rna'][0]['mm1_2'], 0)

        self.assertEqual(
            embedded_data['config']['input']['genome_assembly'],
            'fixture_genome.fa',
        )
        self.assertEqual(
            embedded_data['config']['crispr'],
            {'pam_pattern': 'NAG', 'max_mismatch': 3, 'threads': 7},
        )
        self.assertEqual(
            embedded_data['config']['busco'],
            {'lineage': 'fixture_lineage'},
        )
        self.assertIn('FIXTURE_GENE', html_content)
        self.assertIn('FIXTURE_GUIDE_NAG', html_content)
        self.assertIn('fixture_genome.fa', html_content)
        self.assertIn('<strong>PAM Pattern:</strong> NAG', html_content)
        self.assertIn('<strong>Max Mismatch:</strong> 3', html_content)
        self.assertIn('<strong>Threads:</strong> 7', html_content)

    def test_html_guide_lookup_has_no_gene_name_fallback(self):
        html_content = self.generator.create_pathway_html_template(
            [],
            {'pathway_guide_data': {}, 'pathway_svgs': {}, 'statistics': {}},
            {},
        )

        self.assertIn('guide => guide.gene_key === geneKey', html_content)
        self.assertIn('guide.source_xref_keys.includes(geneKey)', html_content)
        self.assertNotIn('guide.gene_name === geneName', html_content)
        self.assertIn('Pathway Xref Database', html_content)
        self.assertIn('Annotation Source Database', html_content)
        self.assertIn('PW-KFI Gene Key', html_content)
        self.assertIn('Source Xref Keys', html_content)
        self.assertIn('Failure Reason', html_content)
        self.assertIn(
            "guide['MM0'] == null ? '' : guide['MM0']", html_content
        )
        self.assertIn(
            "guide['MM1-2'] == null ? '' : guide['MM1-2']", html_content
        )
        self.assertIn("'MM0',\n                'MM1-2',", html_content)
        self.assertIn('PW-KFI gene coverage', html_content)
        self.assertNotIn('PW-KFI node coverage', html_content)

    @unittest.skipUnless(shutil.which('node'), 'Node.js is not available.')
    def test_guide_csv_exports_mm0_and_mm1_2_values(self):
        html_content = self.generator.create_pathway_html_template(
            [],
            {'pathway_guide_data': {}, 'pathway_svgs': {}, 'statistics': {}},
            {},
        )
        function_marker = 'function generatePathwayGuidesCSVContent'
        function_start = html_content.index(function_marker)
        function_end = html_content.index(
            '\n        }\n    </script>', function_start
        ) + len('\n        }')
        function_source = html_content[function_start:function_end]
        guide_data = {
            'guides': [{
                'pathway_id': 'WP1',
                'MM0': '1',
                'MM1-2': '0',
            }]
        }
        javascript = (
            function_source
            + '\nprocess.stdout.write('
            + 'generatePathwayGuidesCSVContent('
            + json.dumps(guide_data)
            + '));'
        )

        result = subprocess.run(
            [shutil.which('node'), '-e', javascript],
            check=True,
            capture_output=True,
            text=True,
        )
        rows = list(csv.reader(result.stdout.splitlines()))
        mm0_index = rows[0].index('MM0')
        mm1_2_index = rows[0].index('MM1-2')

        self.assertEqual(rows[1][mm0_index], '1')
        self.assertEqual(rows[1][mm1_2_index], '0')


if __name__ == '__main__':
    unittest.main()
