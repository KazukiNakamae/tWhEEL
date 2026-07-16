import json
import shutil
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
            'database': {'force_rebuild': force_rebuild},
            'report': {
                'generate_summary': False,
                'generate_html': True,
            },
            'input': {},
            'crispr': {},
            'busco': {},
        }

    @staticmethod
    def _write_go_kfi(path, go_ids):
        rows = [
            'go_id,total_genes,genes_with_guides,unused,kfi_score,targetable_type'
        ]
        rows.extend(f'{go_id},1,1,,1.0,Fully' for go_id in go_ids)
        path.write_text('\n'.join(rows) + '\n', encoding='utf-8')

    @staticmethod
    def _write_final_guides(path):
        path.write_text(
            '\t'.join([
                'query_id', 'subject_id', 'identity', 'alignment_length',
                'evalue', 'bit_score', 'description', 'gene_name',
                'classification', 'unique_id', 'guide_targets',
            ])
            + '\n'
            + '\t'.join([
                'query1', 'subject1', '99.0', '100', '1e-20', '200',
                'description', 'GENE1', 'Single', 'Single_1',
                'Single_1/10/33/ACGTNGG/1/0/0.987/-',
            ])
            + '\n',
            encoding='utf-8',
        )

    @staticmethod
    def _write_many_final_guides(path, guide_count):
        guide_targets = []
        for index in range(guide_count):
            strand = '+' if index % 2 == 0 else '-'
            guide_targets.append(
                f'Single_{index}/{index}/{index + 23}/'
                f'GUIDE{index:04d}NGG/1/{index % 3}/0.999/{strand}'
            )
        path.write_text(
            '\t'.join([
                'query_id', 'subject_id', 'identity', 'alignment_length',
                'evalue', 'bit_score', 'description', 'gene_name',
                'classification', 'unique_id', 'guide_targets',
            ])
            + '\n'
            + '\t'.join([
                'query1', 'subject1', '99.0', '100', '1e-20', '200',
                'description', 'GENE1', 'Single', 'Single_1',
                ','.join(guide_targets),
            ])
            + '\n',
            encoding='utf-8',
        )

    @staticmethod
    def _results(config, **data_paths):
        return {
            'timestamp': '2026-07-16T00:00:00',
            'statistics': {},
            'config': config,
            'data': {
                key: str(path) for key, path in data_paths.items()
                if path is not None
            },
        }

    @staticmethod
    def _term_info(go_ids):
        return {
            go_id: {
                'name': f'Name for {go_id}',
                'definition': f'Definition for {go_id}',
            }
            for go_id in go_ids
        }

    @staticmethod
    def _go_rows(db_path):
        with sqlite3.connect(db_path) as conn:
            return conn.execute(
                'SELECT go_id, go_name, go_definition '
                'FROM go_terms ORDER BY go_id'
            ).fetchall()

    def _prepare(self, generator, results):
        with mock.patch.object(
            generator,
            '_fetch_go_terms_from_api',
            side_effect=self._term_info,
        ):
            return generator.prepare_web_data(results)

    def test_fresh_database_has_only_the_minimal_cache_schema(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            go_path = temp_dir / 'go_kfi_table.csv'
            self._write_go_kfi(go_path, ['GO:0000001'])
            generator = ReportGenerator(config)

            db_path, _ = self._prepare(
                generator,
                self._results(config, go_kfi_table=go_path),
            )

            with sqlite3.connect(db_path) as conn:
                tables = {
                    row[0] for row in conn.execute(
                        "SELECT name FROM sqlite_master "
                        "WHERE type='table' AND name NOT LIKE 'sqlite_%'"
                    )
                }
                go_columns = [
                    row[1] for row in conn.execute('PRAGMA table_info(go_terms)')
                ]
                metadata = dict(conn.execute(
                    'SELECT key, value FROM report_metadata'
                ))

            self.assertEqual(tables, {'report_metadata', 'go_terms'})
            self.assertEqual(
                go_columns, ['go_id', 'go_name', 'go_definition']
            )
            self.assertEqual(
                set(metadata), {'schema_version', 'input_fingerprint'}
            )
            self.assertEqual(
                metadata['schema_version'],
                str(ReportGenerator.DATABASE_SCHEMA_VERSION),
            )

    def test_matching_go_identifier_set_reuses_database(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            go_path = temp_dir / 'go_kfi_table.csv'
            self._write_go_kfi(go_path, ['GO:0000001'])
            results = self._results(config, go_kfi_table=go_path)
            generator = ReportGenerator(config)

            db_path, _ = self._prepare(generator, results)
            with mock.patch.object(
                generator,
                '_rebuild_database',
                wraps=generator._rebuild_database,
            ) as rebuild_mock:
                reused_path, _ = self._prepare(generator, results)

            self.assertEqual(reused_path, db_path)
            rebuild_mock.assert_not_called()

    def test_standard_run_publishes_one_source_generation(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            go_path = (
                Path(config['output']['intermediate'])
                / 'go_kfi_table.csv'
            )
            self._write_go_kfi(go_path, ['GO:0000001'])
            generator = ReportGenerator(config)

            with mock.patch.object(
                generator.report_source_store,
                'publish_namespace',
                wraps=generator.report_source_store.publish_namespace,
            ) as publish_mock, mock.patch.object(
                generator,
                '_fetch_go_terms_from_api',
                side_effect=self._term_info,
            ):
                generator.run()

            publish_mock.assert_called_once()

    def test_go_identifier_change_rebuilds_without_stale_rows(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            go_path = temp_dir / 'go_kfi_table.csv'
            self._write_go_kfi(go_path, ['GO:0000001'])
            results = self._results(config, go_kfi_table=go_path)
            generator = ReportGenerator(config)

            db_path, _ = self._prepare(generator, results)
            self._write_go_kfi(go_path, ['GO:0000002'])
            rebuilt_generator = ReportGenerator(config)
            rebuilt_path, _ = self._prepare(
                rebuilt_generator,
                self._results(config, go_kfi_table=go_path),
            )

            self.assertEqual(rebuilt_path, db_path)
            self.assertEqual(
                self._go_rows(db_path),
                [(
                    'GO:0000002',
                    'Name for GO:0000002',
                    'Definition for GO:0000002',
                )],
            )

    def test_non_cache_report_source_change_does_not_rebuild_database(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            go_path = temp_dir / 'go_kfi_table.csv'
            guide_path = temp_dir / 'annotation_with_guides_final.tsv'
            self._write_go_kfi(go_path, ['GO:0000001'])
            self._write_final_guides(guide_path)
            results = self._results(
                config,
                go_kfi_table=go_path,
                final_annotation_guides=guide_path,
            )
            generator = ReportGenerator(config)
            self._prepare(generator, results)
            guide_path.write_text(
                guide_path.read_text(encoding='utf-8') + '\n',
                encoding='utf-8',
            )
            next_generator = ReportGenerator(config)
            next_results = self._results(
                config,
                go_kfi_table=go_path,
                final_annotation_guides=guide_path,
            )

            with mock.patch.object(
                next_generator,
                '_rebuild_database',
                wraps=next_generator._rebuild_database,
            ) as rebuild_mock:
                self._prepare(next_generator, next_results)

            rebuild_mock.assert_not_called()

    def test_force_rebuild_replaces_cache_instead_of_appending(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            go_path = temp_dir / 'go_kfi_table.csv'
            self._write_go_kfi(go_path, ['GO:0000001'])
            results = self._results(config, go_kfi_table=go_path)
            generator = ReportGenerator(config)

            db_path, _ = self._prepare(generator, results)
            generator.config['database']['force_rebuild'] = True
            self._prepare(generator, results)

            self.assertEqual(len(self._go_rows(db_path)), 1)

    def test_empty_cache_is_reused_when_go_source_is_absent(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            results = self._results(config)
            generator = ReportGenerator(config)

            db_path, _ = self._prepare(generator, results)
            with mock.patch.object(
                generator,
                '_rebuild_database',
                wraps=generator._rebuild_database,
            ) as rebuild_mock:
                self._prepare(generator, results)

            rebuild_mock.assert_not_called()
            self.assertEqual(self._go_rows(db_path), [])

    def test_failed_force_rebuild_preserves_published_cache(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            go_path = temp_dir / 'go_kfi_table.csv'
            self._write_go_kfi(go_path, ['GO:0000001'])
            results = self._results(config, go_kfi_table=go_path)
            generator = ReportGenerator(config)

            db_path, _ = self._prepare(generator, results)
            original_rows = self._go_rows(db_path)
            generator.config['database']['force_rebuild'] = True
            with mock.patch.object(
                generator,
                '_insert_go_data',
                side_effect=RuntimeError('insertion failed'),
            ):
                with self.assertRaisesRegex(RuntimeError, 'insertion failed'):
                    generator.prepare_web_data(results)

            self.assertEqual(self._go_rows(db_path), original_rows)
            self.assertEqual(
                list(Path(generator.web_data_dir).glob('.twheel_data.*.tmp')),
                [],
            )

    def test_legacy_wide_database_is_rebuilt_and_reuses_go_labels(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            go_path = temp_dir / 'go_kfi_table.csv'
            self._write_go_kfi(go_path, ['GO:0000001'])
            generator = ReportGenerator(config)
            Path(generator.web_data_dir).mkdir(parents=True, exist_ok=True)
            db_path = Path(generator.web_data_dir) / 'twheel_data.db'
            with sqlite3.connect(db_path) as conn:
                conn.execute('CREATE TABLE genes (id INTEGER)')
                conn.execute('''
                    CREATE TABLE go_terms (
                        id INTEGER PRIMARY KEY,
                        go_id TEXT,
                        go_name TEXT,
                        go_definition TEXT,
                        total_genes INTEGER
                    )
                ''')
                conn.execute(
                    'INSERT INTO go_terms '
                    '(go_id, go_name, go_definition, total_genes) '
                    "VALUES ('GO:0000001', 'Legacy name', 'Legacy definition', 1)"
                )

            with mock.patch.object(
                generator, '_fetch_go_terms_from_api'
            ) as fetch_mock:
                rebuilt_path, _ = generator.prepare_web_data(
                    self._results(config, go_kfi_table=go_path)
                )

            fetch_mock.assert_not_called()
            self.assertEqual(
                self._go_rows(rebuilt_path),
                [('GO:0000001', 'Legacy name', 'Legacy definition')],
            )
            with sqlite3.connect(rebuilt_path) as conn:
                tables = {
                    row[0] for row in conn.execute(
                        "SELECT name FROM sqlite_master "
                        "WHERE type='table' AND name NOT LIKE 'sqlite_%'"
                    )
                }
            self.assertEqual(tables, {'report_metadata', 'go_terms'})

    def test_report_regenerates_after_intermediate_directory_is_removed(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            intermediate = Path(config['output']['intermediate'])
            go_path = intermediate / 'go_kfi_table.csv'
            guide_path = intermediate / 'annotation_with_guides_final.tsv'
            reverse_path = intermediate / 'go_reverse_lookup.tsv'
            parent_path = intermediate / 'go_parent_counts.csv'
            image_path = intermediate / 'go_parent_counts.png'
            treemap_path = intermediate / 'treemap_go_parent_counts.png'
            (intermediate / 'transcriptome.fa').write_text(
                '>transcript_1\nATGCATGC\n', encoding='utf-8'
            )
            (intermediate / 'query.pep').write_text(
                '>protein_1\nMPEPTIDE\n', encoding='utf-8'
            )
            (intermediate / 'annotation_classified.txt').write_text(
                '\t'.join([
                    'query_id', 'subject_id', 'identity', 'alignment_length',
                    'evalue', 'bit_score', 'description', 'gene_name',
                    'classification', 'unique_id',
                ])
                + '\n'
                + '\t'.join([
                    'query1', 'subject1', '99.0', '100', '1e-20', '200',
                    'description', 'GENE1', 'Single', 'Single_1',
                ])
                + '\n',
                encoding='utf-8',
            )
            (intermediate / 'multi_genes_targets_specific_clean.bed').write_text(
                'Multi_1\t1\t24\n', encoding='utf-8'
            )
            (intermediate / 'single_genes_targets_specific_clean.bed').write_text(
                'Single_1\t10\t33\n', encoding='utf-8'
            )
            (intermediate / 'go_kfi_stat.csv').write_text(
                'F:100%[S:100%,M:0%],P:0%,N:0%,n:1\n'
                '1 Fully targetable GOs (F)\n',
                encoding='utf-8',
            )
            self._write_go_kfi(go_path, ['GO:0000001'])
            self._write_final_guides(guide_path)
            reverse_path.write_text(
                'go_term\tgene_count\tgenes\n'
                'GO:0000001\t1\tSingle_1\n',
                encoding='utf-8',
            )
            parent_path.write_text(
                'parent_go,description,count\n'
                'GO:0008150,biological_process,1\n',
                encoding='utf-8',
            )
            image_path.write_bytes(b'not-a-real-png')
            treemap_path.write_bytes(b'not-a-real-treemap-png')
            generator = ReportGenerator(config)
            with mock.patch.object(
                generator,
                '_fetch_go_terms_from_api',
                side_effect=self._term_info,
            ):
                first_results = generator.run()

            staged_guide = Path(
                first_results['report_data']['final_annotation_guides']
            )
            first_statistics = first_results['statistics']
            first_manifest = json.loads(
                Path(first_results['report_sources_manifest']).read_text(
                    encoding='utf-8'
                )
            )
            self.assertTrue(staged_guide.is_file())
            self.assertEqual(
                set(first_results['report_data']),
                set(ReportGenerator.REPORT_SOURCE_FILES),
            )
            self.assertEqual(
                {
                    key for key in first_manifest['sources']
                    if key.startswith('normal.')
                },
                {
                    f'normal.{key}'
                    for key in ReportGenerator.REPORT_SOURCE_FILES
                },
            )
            shutil.rmtree(intermediate)

            regenerated = ReportGenerator(config).run()
            html = Path(regenerated['html_path']).read_text(encoding='utf-8')
            marker = '<script type="application/json" id="twheel-data">'
            embedded_json = html.split(marker, 1)[1].split('</script>', 1)[0]
            embedded_data = json.loads(embedded_json)

            self.assertTrue(staged_guide.is_file())
            self.assertEqual(len(embedded_data['gene_annotation']), 1)
            self.assertEqual(len(embedded_data['guide_rna']), 1)
            self.assertEqual(len(embedded_data['go_analysis']), 1)
            self.assertEqual(len(embedded_data['go_parent_counts']), 1)
            self.assertEqual(regenerated['statistics'], first_statistics)
            self.assertEqual(
                embedded_data['statistics'], first_statistics
            )
            self.assertEqual(
                regenerated['statistics']['gene_prediction']['transcript_count'],
                1,
            )
            self.assertEqual(
                regenerated['statistics']['gene_prediction']['protein_count'],
                1,
            )
            self.assertEqual(
                regenerated['statistics']['annotation']['total_genes'], 1
            )
            self.assertEqual(
                regenerated['statistics']['guide_design']['multi_targets'], 1
            )
            self.assertEqual(
                regenerated['statistics']['guide_design']['single_targets'], 1
            )
            self.assertEqual(
                regenerated['statistics']['go_analysis']['total_gos'], 1
            )
            self.assertEqual(
                (Path(generator.output_dir) / 'static' / 'img'
                 / 'go_parent_counts.png').read_bytes(),
                b'not-a-real-png',
            )
            self.assertEqual(
                (Path(generator.output_dir) / 'static' / 'img'
                 / 'treemap_go_parent_counts.png').read_bytes(),
                b'not-a-real-treemap-png',
            )
            self.assertTrue(
                Path(regenerated['report_sources_manifest']).is_file()
            )
            self.assertTrue(
                Path(regenerated['report_generation_results_path']).is_file()
            )
            self.assertFalse(intermediate.exists())

    def test_html_embeds_all_guides_beyond_first_thousand(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            final_guides_path = temp_dir / 'many_final_guides.tsv'
            self._write_many_final_guides(final_guides_path, 1005)
            results = self._results(
                config, final_annotation_guides=final_guides_path
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
                embedded_data['guide_rna'][-1]['reference'], 'Single_1004'
            )

            generator._create_js_files()
            javascript = (
                Path(generator.output_dir) / 'static' / 'js' / 'main.js'
            ).read_text(encoding='utf-8')
            self.assertNotIn('const rows = guideData.map', javascript)
            self.assertIn('this.setupGuideTableFilters(guideData);', javascript)

    def test_new_generation_prunes_missing_optional_images(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            intermediate = Path(config['output']['intermediate'])
            go_path = intermediate / 'go_kfi_table.csv'
            image_path = intermediate / 'go_parent_counts.png'
            treemap_path = intermediate / 'treemap_go_parent_counts.png'
            self._write_go_kfi(go_path, ['GO:0000001'])
            image_path.write_bytes(b'old-chart')
            treemap_path.write_bytes(b'old-treemap')

            generator = ReportGenerator(config)
            with mock.patch.object(
                generator,
                '_fetch_go_terms_from_api',
                side_effect=self._term_info,
            ):
                generator.run()

            image_path.unlink()
            treemap_path.unlink()
            regenerated = ReportGenerator(config).run()
            manifest = json.loads(
                Path(regenerated['report_sources_manifest']).read_text(
                    encoding='utf-8'
                )
            )

            self.assertNotIn(
                'go_parent_counts_image', regenerated['report_data']
            )
            self.assertNotIn(
                'go_parent_counts_treemap_image',
                regenerated['report_data'],
            )
            self.assertNotIn(
                'normal.go_parent_counts_image', manifest['sources']
            )
            self.assertNotIn(
                'normal.go_parent_counts_treemap_image',
                manifest['sources'],
            )
            self.assertFalse(
                (Path(generator.output_dir) / 'static' / 'img'
                 / 'go_parent_counts.png').exists()
            )
            self.assertFalse(
                (Path(generator.output_dir) / 'static' / 'img'
                 / 'treemap_go_parent_counts.png').exists()
            )

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
