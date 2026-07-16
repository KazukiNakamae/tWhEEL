import json
import sqlite3
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

try:
    import requests  # noqa: F401
except ModuleNotFoundError:
    sys.modules['requests'] = types.ModuleType('requests')

from src.modules.report_generator import ReportGenerator


class ReportDatabaseTests(unittest.TestCase):
    @staticmethod
    def _config(temp_dir, force_rebuild=False):
        output_dir = temp_dir / 'output'
        intermediate_dir = output_dir / 'intermediate'
        intermediate_dir.mkdir(parents=True, exist_ok=True)
        return {
            'output': {
                'base_dir': str(output_dir),
                'intermediate': str(intermediate_dir),
            },
            'database': {
                'force_rebuild': force_rebuild,
            },
            'input': {},
            'crispr': {},
            'busco': {},
        }

    @staticmethod
    def _write_annotation(annotation_path, query_id):
        annotation_path.write_text(
            '\t'.join([
                'query_id',
                'subject_id',
                'identity',
                'alignment_length',
                'evalue',
                'bit_score',
                'description',
                'gene_name',
                'classification',
                'unique_id',
            ])
            + '\n'
            + '\t'.join([
                query_id,
                'subject1',
                '99.0',
                '100',
                '1e-20',
                '200',
                'description',
                'GENE1',
                'Single',
                'Single_1',
            ])
            + '\n',
            encoding='utf-8',
        )

    @staticmethod
    def _write_final_guides(final_guides_path):
        final_guides_path.write_text(
            '\t'.join([
                'query_id',
                'subject_id',
                'identity',
                'alignment_length',
                'evalue',
                'bit_score',
                'description',
                'gene_name',
                'classification',
                'unique_id',
                'guide_targets',
            ])
            + '\n'
            + '\t'.join([
                'query1',
                'subject1',
                '99.0',
                '100',
                '1e-20',
                '200',
                'description',
                'GENE1',
                'Single',
                'Single_1',
                ','.join([
                    'Single_1/10/33/ACGTNGG/1/0/0.987/-',
                    'Single_1/20/43/TGCANGG/1/2',
                    'Single_1/30/53/GGCANGG/1/0/0.999/invalid',
                ]),
            ])
            + '\n',
            encoding='utf-8',
        )

    @staticmethod
    def _write_many_final_guides(final_guides_path, guide_count):
        guide_targets = []
        for index in range(guide_count):
            strand = '+' if index % 2 == 0 else '-'
            guide_targets.append(
                f'Single_{index}/{index}/{index + 23}/'
                f'GUIDE{index:04d}NGG/1/{index % 3}/0.999/{strand}'
            )

        final_guides_path.write_text(
            '\t'.join([
                'query_id',
                'subject_id',
                'identity',
                'alignment_length',
                'evalue',
                'bit_score',
                'description',
                'gene_name',
                'classification',
                'unique_id',
                'guide_targets',
            ])
            + '\n'
            + '\t'.join([
                'query1',
                'subject1',
                '99.0',
                '100',
                '1e-20',
                '200',
                'description',
                'GENE1',
                'Single',
                'Single_1',
                ','.join(guide_targets),
            ])
            + '\n',
            encoding='utf-8',
        )

    @staticmethod
    def _results(config, annotation_path=None, final_guides_path=None):
        data = {}
        if annotation_path is not None:
            data['annotation_table'] = str(annotation_path)
        if final_guides_path is not None:
            data['final_annotation_guides'] = str(final_guides_path)
        return {
            'timestamp': '2026-07-15T00:00:00',
            'statistics': {},
            'config': config,
            'data': data,
        }

    @staticmethod
    def _gene_query_ids(db_path):
        with sqlite3.connect(db_path) as conn:
            return [
                row[0]
                for row in conn.execute(
                    'SELECT query_id FROM genes ORDER BY id'
                ).fetchall()
            ]

    def test_matching_source_database_is_reused(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            annotation_path = temp_dir / 'annotation.tsv'
            self._write_annotation(annotation_path, 'query1')
            results = self._results(config, annotation_path)
            generator = ReportGenerator(config)

            db_path, _ = generator.prepare_web_data(results)
            with mock.patch.object(
                generator,
                '_rebuild_database',
                wraps=generator._rebuild_database,
            ) as rebuild_mock:
                reused_path, _ = generator.prepare_web_data(results)

            self.assertEqual(reused_path, db_path)
            rebuild_mock.assert_not_called()
            self.assertEqual(self._gene_query_ids(db_path), ['query1'])

    def test_source_content_change_rebuilds_without_stale_rows(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            annotation_path = temp_dir / 'annotation.tsv'
            self._write_annotation(annotation_path, 'query1')
            results = self._results(config, annotation_path)
            generator = ReportGenerator(config)

            db_path, _ = generator.prepare_web_data(results)
            self._write_annotation(annotation_path, 'query2')
            rebuilt_path, _ = generator.prepare_web_data(results)

            self.assertEqual(rebuilt_path, db_path)
            self.assertEqual(self._gene_query_ids(db_path), ['query2'])

    def test_force_rebuild_replaces_database_instead_of_appending(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            annotation_path = temp_dir / 'annotation.tsv'
            self._write_annotation(annotation_path, 'query1')
            results = self._results(config, annotation_path)
            generator = ReportGenerator(config)

            db_path, _ = generator.prepare_web_data(results)
            generator.config['database']['force_rebuild'] = True
            generator.prepare_web_data(results)

            self.assertEqual(self._gene_query_ids(db_path), ['query1'])

    def test_empty_database_is_reused_when_empty_sources_match(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            results = self._results(config)
            generator = ReportGenerator(config)

            db_path, _ = generator.prepare_web_data(results)
            with mock.patch.object(
                generator,
                '_rebuild_database',
                wraps=generator._rebuild_database,
            ) as rebuild_mock:
                generator.prepare_web_data(results)

            rebuild_mock.assert_not_called()
            self.assertEqual(self._gene_query_ids(db_path), [])

    def test_failed_force_rebuild_preserves_published_database(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            annotation_path = temp_dir / 'annotation.tsv'
            self._write_annotation(annotation_path, 'query1')
            results = self._results(config, annotation_path)
            generator = ReportGenerator(config)

            db_path, _ = generator.prepare_web_data(results)
            generator.config['database']['force_rebuild'] = True
            with mock.patch.object(
                generator,
                '_insert_gene_data',
                side_effect=RuntimeError('insertion failed'),
            ):
                with self.assertRaisesRegex(RuntimeError, 'insertion failed'):
                    generator.prepare_web_data(results)

            self.assertEqual(self._gene_query_ids(db_path), ['query1'])
            self.assertEqual(
                list(Path(generator.web_data_dir).glob('.twheel_data.*.tmp')),
                [],
            )

    def test_guide_strand_uses_final_field_instead_of_mm1_2(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            annotation_path = temp_dir / 'annotation.tsv'
            final_guides_path = temp_dir / 'final_guides.tsv'
            self._write_annotation(annotation_path, 'query1')
            self._write_final_guides(final_guides_path)
            results = self._results(
                config,
                annotation_path,
                final_guides_path,
            )
            generator = ReportGenerator(config)

            db_path, _ = generator.prepare_web_data(results)

            with sqlite3.connect(db_path) as conn:
                guides = conn.execute(
                    "SELECT chromosome, start, end, sequence, strand "
                    "FROM guide_rnas ORDER BY id"
                ).fetchall()
                guide_count = conn.execute(
                    "SELECT guide_count FROM genes WHERE query_id = 'query1'"
                ).fetchone()[0]

            self.assertEqual(
                guides,
                [('Single_1', 10, 33, 'ACGTNGG', '-')],
            )
            self.assertEqual(guide_count, 1)

    def test_html_embeds_all_guides_beyond_first_thousand(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            final_guides_path = temp_dir / 'many_final_guides.tsv'
            self._write_many_final_guides(final_guides_path, 1005)
            results = self._results(
                config,
                final_guides_path=final_guides_path,
            )
            generator = ReportGenerator(config)

            guide_data = generator._load_guide_rna_data_for_js(results)
            html = generator._generate_html_content(results)
            marker = '<script type="application/json" id="twheel-data">'
            embedded_json = html.split(marker, 1)[1].split('</script>', 1)[0]
            embedded_data = json.loads(embedded_json)

            self.assertEqual(len(guide_data), 1005)
            self.assertEqual(guide_data[-1]['reference'], 'Single_1004')
            self.assertEqual(guide_data[-1]['target_sequence'], 'GUIDE1004NGG')
            self.assertEqual(len(embedded_data['guide_rna']), 1005)
            self.assertEqual(
                embedded_data['guide_rna'][-1]['reference'],
                'Single_1004',
            )

            generator._create_js_files()
            javascript = (
                Path(generator.output_dir) / 'static' / 'js' / 'main.js'
            ).read_text(encoding='utf-8')
            self.assertNotIn('const rows = guideData.map', javascript)
            self.assertIn('this.setupGuideTableFilters(guideData);', javascript)

    def test_go_reverse_lookup_accepts_normalized_and_legacy_entries(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            reverse_lookup_path = temp_dir / 'go_reverse_lookup.tsv'
            reverse_lookup_path.write_text(
                'go_term\tgene_count\tgenes\n'
                'GO:0000001\t2\tSingle_1;Single_2\n'
                'GO:0000002\t1\tSingle_3|description|guide_target\n',
                encoding='utf-8',
            )
            generator = ReportGenerator(self._config(temp_dir))

            go_genes = generator._parse_go_annotated_genes({
                'data': {'go_reverse_lookup': str(reverse_lookup_path)},
            })

        self.assertEqual(
            go_genes,
            {
                'GO:0000001': ['Single_1', 'Single_2'],
                'GO:0000002': ['Single_3'],
            },
        )


if __name__ == '__main__':
    unittest.main()
