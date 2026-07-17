import json
import re
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
    def _browser_report_data(output_dir, report_name='report.html'):
        """Load the external core payload and reconstruct deferred guides."""
        output_dir = Path(output_dir)
        html = (output_dir / report_name).read_text(encoding='utf-8')
        match = re.search(
            r'<script src="([^"]*report_payloads/[^"]+/core\.js)"></script>',
            html,
        )
        if match is None:
            raise AssertionError('Report HTML does not reference a payload')
        core_path = output_dir / match.group(1)
        core_source = core_path.read_text(encoding='utf-8')
        marker = 'var payload = '
        payload, _ = json.JSONDecoder().raw_decode(
            core_source[core_source.index(marker) + len(marker):]
        )
        data = dict(payload['data'])
        guides = []
        for chunk in payload['guides']['chunks']:
            chunk_source = (core_path.parent / chunk['src']).read_text(
                encoding='utf-8'
            )
            assignment = chunk_source.index('] = ') + len('] = ')
            rows, _ = json.JSONDecoder().raw_decode(
                chunk_source[assignment:]
            )
            guides.extend(rows)
        if len(guides) != payload['guides']['row_count']:
            raise AssertionError('Guide payload row count mismatch')
        data['guide_rna'] = guides
        return html, payload, data

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

    def test_report_scope_warning_distinguishes_genome_reference_mode(self):
        scenarios = (
            (
                '',
                'transcriptome_only_preliminary',
                'Preliminary feasibility screening only',
                (
                    'splice-junction',
                    'does not establish a contiguous genomic target',
                    'mismatch limit fixed at 2',
                ),
            ),
            (
                '/analysis/genome.fa',
                'genome_off_target_reference',
                'Genome-backed target-locus validation required',
                (
                    'up to 2 mismatches',
                    'MM0=0',
                    'MM0=1',
                    'MM1-2=0',
                    'intended gene locus',
                    'GFF/GTF',
                ),
            ),
        )

        for genome_path, mode, title, expected_phrases in scenarios:
            with self.subTest(mode=mode):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    config = self._config(Path(temp_dir_name))
                    if genome_path:
                        config['input']['genome_assembly'] = genome_path
                    generator = ReportGenerator(config)
                    results = self._results(config)

                    html = generator._generate_html_content(
                        results,
                        payload_script_src=(
                            'web_data/report_payloads/test/core.js'
                        ),
                    )
                    generator._create_css_files()
                    css = (
                        Path(generator.output_dir)
                        / 'static'
                        / 'css'
                        / 'main.css'
                    ).read_text(encoding='utf-8')
                    summary_path = generator.generate_text_summary(results)
                    summary = Path(summary_path).read_text(encoding='utf-8')

                marker = f'data-guide-design-scope="{mode}"'
                self.assertEqual(html.count(marker), 1)
                self.assertIn(title, html)
                for expected_phrase in expected_phrases:
                    self.assertIn(expected_phrase, html)
                self.assertLess(
                    html.index(marker),
                    html.index('<!-- General Statistics Section -->'),
                )
                self.assertIn(title, summary)
                for expected_phrase in expected_phrases:
                    self.assertIn(expected_phrase, summary)
                self.assertIn(
                    '<strong>Maximum Mismatch (fixed by tWhEEL):</strong> 2',
                    html,
                )
                self.assertIn(
                    'Maximum mismatch (fixed by tWhEEL): 2', summary
                )
                if genome_path:
                    self.assertNotIn(
                        'otherwise lack a contiguous genomic target', html
                    )
                    self.assertNotIn(
                        'otherwise lack a contiguous genomic target', summary
                    )
                self.assertIn('.alert-critical', css)
                self.assertIn('border: 2px solid #b02a37', css)

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
                run_results = generator.run()

            publish_mock.assert_called_once()
            summary = json.loads(
                Path(run_results['json_path']).read_text(encoding='utf-8')
            )
            browser_pointer = summary['browser_payload']
            self.assertEqual(browser_pointer['schema_version'], 2)
            self.assertEqual(browser_pointer['namespace'], 'normal')
            self.assertTrue(
                (Path(generator.output_dir) / browser_pointer['manifest'])
                .is_file()
            )
            self.assertTrue(
                (
                    Path(generator.web_data_dir)
                    / 'report_payloads'
                    / '.committed'
                    / Path(run_results['browser_payload']['core_path']).parent.name
                ).is_file()
            )
            report_html = Path(run_results['html_path']).read_text(
                encoding='utf-8'
            )
            self.assertNotIn('type="application/json"', report_html)

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

    def test_normalized_report_sources_replace_legacy_guide_table_snapshot(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            legacy_path = temp_dir / 'annotation_with_guides_final.tsv'
            genes_path = temp_dir / 'normalized_genes.tsv'
            guides_path = temp_dir / 'normalized_guides.tsv'
            manifest_path = temp_dir / 'normalized_gene_guide_manifest.json'
            self._write_final_guides(legacy_path)
            genes_path.write_text('normalized genes\n', encoding='utf-8')
            guides_path.write_text('normalized guides\n', encoding='utf-8')
            manifest_path.write_text('{}\n', encoding='utf-8')
            results = self._results(
                config,
                final_annotation_guides=legacy_path,
                normalized_genes=genes_path,
                normalized_guides=guides_path,
                normalized_gene_guide_manifest=manifest_path,
            )

            generator = ReportGenerator(config)
            generator._stage_report_sources(results)

            self.assertNotIn(
                'final_annotation_guides', results['report_data']
            )
            self.assertEqual(
                set(ReportGenerator.NORMALIZED_GENE_GUIDE_SOURCE_FILES),
                set(results['report_data']),
            )
            for path in results['report_data'].values():
                self.assertTrue(Path(path).is_file())

    def test_partial_normalized_sources_keep_legacy_snapshot(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            legacy_path = temp_dir / 'annotation_with_guides_final.tsv'
            genes_path = temp_dir / 'normalized_genes.tsv'
            self._write_final_guides(legacy_path)
            genes_path.write_text('partial normalized data\n', encoding='utf-8')
            results = self._results(
                config,
                final_annotation_guides=legacy_path,
                normalized_genes=genes_path,
            )

            ReportGenerator(config)._stage_report_sources(results)

            self.assertIn('final_annotation_guides', results['report_data'])
            self.assertNotIn('normalized_genes', results['report_data'])

    def test_normalized_gene_and_guide_data_are_loaded_only_once(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            genes_path = temp_dir / 'normalized_genes.tsv'
            guides_path = temp_dir / 'normalized_guides.tsv'
            manifest_path = temp_dir / 'normalized_gene_guide_manifest.json'
            for path in (genes_path, guides_path, manifest_path):
                path.write_text('cache-key input\n', encoding='utf-8')

            genes = [
                    {
                        'unique_id': 'Multi_1',
                        'classification': 'Multi',
                        'subject_id': 'subject_multi',
                        'gene_name': 'MULTI',
                        'description': 'multi description',
                        'representative_query_id': 'contig_a',
                        'query_ids': ['contig_a', 'contig_b'],
                        'annotation_row_count': 2,
                        'guide_count': 2,
                        'has_guides': True,
                    },
                    {
                        'unique_id': 'Single_1',
                        'classification': 'Single',
                        'subject_id': 'subject_single',
                        'gene_name': 'SINGLE',
                        'description': 'single description',
                        'representative_query_id': 'contig_c',
                        'query_ids': ['contig_c'],
                        'annotation_row_count': 1,
                        'guide_count': 1,
                        'has_guides': True,
                    },
                ]
            guides = [
                    {
                        'unique_id': 'Multi_1',
                        'reference': 'Multi_1',
                        'start': 10,
                        'end': 33,
                        'guide_sequence': 'MULTI_A',
                        'mm0': 1,
                        'mm1_2': 0,
                        'off_target_score': 1.0,
                        'strand': '+',
                    },
                    {
                        'unique_id': 'Multi_1',
                        'reference': 'Multi_1',
                        'start': 20,
                        'end': 43,
                        'guide_sequence': 'MULTI_B',
                        'mm0': 1,
                        'mm1_2': 0,
                        'off_target_score': None,
                        'strand': '',
                    },
                    {
                        'unique_id': 'Single_1',
                        'reference': 'Single_1',
                        'start': 30,
                        'end': 53,
                        'guide_sequence': 'SINGLE',
                        'mm0': 0,
                        'mm1_2': 0,
                        'off_target_score': 0.5,
                        'strand': '-',
                    },
                ]
            loader = mock.Mock(return_value={
                'genes': genes,
                'paths': {
                    'normalized_genes': str(genes_path),
                    'normalized_guides': str(guides_path),
                    'normalized_gene_guide_manifest': str(manifest_path),
                },
            })
            guide_iterator = mock.Mock(
                side_effect=lambda *args, **kwargs: iter(guides)
            )
            normalization_module = types.ModuleType(
                'src.modules.gene_guide_normalization'
            )
            normalization_module.load_normalized_gene_data = loader
            normalization_module.iter_normalized_guides = guide_iterator
            results = self._results(
                config,
                normalized_genes=genes_path,
                normalized_guides=guides_path,
                normalized_gene_guide_manifest=manifest_path,
            )
            generator = ReportGenerator(config)

            with mock.patch.dict(sys.modules, {
                'src.modules.gene_guide_normalization': normalization_module,
            }):
                statistics = generator._collect_statistics(results)
                self.assertIsNone(
                    generator._normalized_guide_projection_cache
                )
                gene_data = generator._load_gene_annotation_data(results)
                guide_data = generator._load_guide_rna_data_for_js(results)

            loader.assert_called_once_with(
                str(genes_path), str(guides_path), str(manifest_path)
            )
            self.assertEqual(guide_iterator.call_count, 2)
            guide_iterator.assert_has_calls([
                mock.call(
                    str(guides_path), str(manifest_path), verify=False
                ),
                mock.call(
                    str(guides_path), str(manifest_path), verify=False
                ),
            ])
            self.assertEqual(statistics['annotation'], {
                'multi_genes': 1,
                'single_genes': 1,
                'total_genes': 2,
            })
            self.assertEqual(statistics['guide_design'], {
                'total_guides': 3,
                'site_specific_guides': 2,
                'site_specific_multi_target': 2,
                'site_specific_single_target': 0,
            })
            self.assertEqual(gene_data[0]['query_id'], 'contig_a; contig_b')
            self.assertEqual(gene_data[0]['target_count'], 2)
            self.assertEqual(guide_data[1]['off_target_score'], 0.0)
            self.assertEqual(guide_data[2]['target_sequence'], 'SINGLE')

    def test_current_normalized_pointer_resolves_exact_generation_sources(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            intermediate = Path(config['output']['intermediate'])
            pointer_path = (
                intermediate
                / ReportGenerator.NORMALIZED_GENE_GUIDE_CURRENT_FILENAME
            )
            pointer_path.write_text('{}\n', encoding='utf-8')
            sources = {}
            for key, filename in (
                ReportGenerator.NORMALIZED_GENE_GUIDE_SOURCE_FILES.items()
            ):
                path = temp_dir / 'generation' / filename
                path.parent.mkdir(exist_ok=True)
                path.write_text(f'{key}\n', encoding='utf-8')
                sources[key] = str(path)

            resolver = mock.Mock(return_value=sources)
            normalization_module = types.ModuleType(
                'src.modules.gene_guide_normalization'
            )
            normalization_module.resolve_current_normalized_gene_guide_paths = (
                resolver
            )
            generator = ReportGenerator(config)

            with mock.patch.dict(sys.modules, {
                'src.modules.gene_guide_normalization': normalization_module,
            }):
                resolved = generator._current_normalized_gene_guide_sources()

            resolver.assert_called_once_with(str(intermediate))
            self.assertEqual(resolved, sources)

    def test_invalid_current_normalized_pointer_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            intermediate = Path(config['output']['intermediate'])
            (
                intermediate
                / ReportGenerator.NORMALIZED_GENE_GUIDE_CURRENT_FILENAME
            ).write_text('{}\n', encoding='utf-8')
            normalization_module = types.ModuleType(
                'src.modules.gene_guide_normalization'
            )
            normalization_module.resolve_current_normalized_gene_guide_paths = (
                mock.Mock(side_effect=ValueError('invalid current pointer'))
            )

            with mock.patch.dict(sys.modules, {
                'src.modules.gene_guide_normalization': normalization_module,
            }):
                with self.assertRaisesRegex(
                    ValueError, 'invalid current pointer'
                ):
                    ReportGenerator(
                        config
                    )._current_normalized_gene_guide_sources()

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
            expected_legacy_sources = (
                set(ReportGenerator.REPORT_SOURCE_FILES)
                - set(ReportGenerator.NORMALIZED_GENE_GUIDE_SOURCE_FILES)
            )
            self.assertEqual(
                set(first_results['report_data']),
                expected_legacy_sources,
            )
            self.assertEqual(
                {
                    key for key in first_manifest['sources']
                    if key.startswith('normal.')
                },
                {
                    f'normal.{key}'
                    for key in expected_legacy_sources
                },
            )
            shutil.rmtree(intermediate)

            regenerated = ReportGenerator(config).run()
            html, browser_payload, embedded_data = (
                self._browser_report_data(generator.output_dir)
            )

            self.assertTrue(staged_guide.is_file())
            self.assertEqual(len(embedded_data['gene_annotation']), 1)
            self.assertEqual(len(embedded_data['guide_rna']), 1)
            self.assertEqual(len(embedded_data['go_analysis']), 1)
            self.assertEqual(len(embedded_data['go_parent_counts']), 1)
            self.assertEqual(regenerated['statistics'], first_statistics)
            self.assertEqual(
                embedded_data['statistics'], first_statistics
            )
            self.assertEqual(browser_payload['schema_version'], 2)
            self.assertNotIn('twheel-data', html)
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

    def test_normalized_snapshot_regenerates_without_legacy_guide_table(self):
        from src.modules.gene_guide_normalization import (
            publish_normalized_gene_guide_generation,
        )

        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            intermediate = Path(config['output']['intermediate'])
            annotation_path = intermediate / 'annotation_classified.txt'
            multi_bed_path = intermediate / 'multi.bed'
            single_bed_path = intermediate / 'single.bed'
            annotation_path.write_text(
                '\t'.join([
                    'query_id', 'subject_id', 'identity', 'length',
                    'evalue', 'bit_score', 'description', 'gene_name',
                    'classification', 'unique_id',
                ])
                + '\n'
                + '\t'.join([
                    'contig_a', 'subject_multi', '99', '100', '1e-30',
                    '200', 'multi description', 'MULTI', 'Multi', 'Multi_1',
                ])
                + '\n'
                + '\t'.join([
                    'contig_b', 'subject_multi', '98', '90', '1e-20',
                    '180', 'multi description', 'MULTI', 'Multi', 'Multi_1',
                ])
                + '\n'
                + '\t'.join([
                    'contig_c', 'subject_single', '97', '80', '1e-10',
                    '160', 'single description', 'SINGLE', 'Single',
                    'Single_1',
                ])
                + '\n',
                encoding='utf-8',
            )
            multi_bed_path.write_text(
                'Multi_1_aligned_consensus\t10\t33\tMULTI_GUIDE/1/0\t1.0\t+\n',
                encoding='utf-8',
            )
            single_bed_path.write_text(
                'Single_1\t20\t43\tSINGLE_GUIDE/1/0\t0.9\t-\n',
                encoding='utf-8',
            )
            publish_normalized_gene_guide_generation(
                annotation_path=annotation_path,
                output_dir=intermediate,
                multi_bed_path=multi_bed_path,
                single_bed_path=single_bed_path,
            )

            first_generator = ReportGenerator(config)
            first_results = first_generator.run()
            self.assertIsNone(
                first_generator._normalized_guide_projection_cache
            )
            self.assertNotIn(
                'final_annotation_guides', first_results['report_data']
            )
            self.assertTrue(
                set(ReportGenerator.NORMALIZED_GENE_GUIDE_SOURCE_FILES)
                <= set(first_results['report_data'])
            )
            shutil.rmtree(intermediate)

            regenerated_generator = ReportGenerator(config)
            regenerated = regenerated_generator.run()
            self.assertIsNone(
                regenerated_generator._normalized_guide_projection_cache
            )
            html, browser_payload, embedded_data = (
                self._browser_report_data(regenerated_generator.output_dir)
            )

            self.assertNotIn(
                'final_annotation_guides', regenerated['report_data']
            )
            self.assertEqual(len(embedded_data['gene_annotation']), 2)
            self.assertEqual(len(embedded_data['guide_rna']), 2)
            self.assertEqual(
                embedded_data['statistics']['annotation'],
                {'multi_genes': 1, 'single_genes': 1, 'total_genes': 2},
            )
            self.assertEqual(
                embedded_data['statistics']['guide_design'][
                    'site_specific_guides'
                ],
                2,
            )
            self.assertEqual(browser_payload['guides']['row_count'], 2)
            self.assertNotIn('twheel-data', html)

    def test_payload_chunks_all_guides_beyond_first_thousand(self):
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
            generator._create_js_files()
            with generator.report_payload_store.publication_lock():
                published = generator.publish_browser_payload(results)
                html = generator._generate_html_content(
                    results,
                    payload_script_src=published['core_relative_path'],
                )
                generator._atomic_write_text(
                    Path(generator.output_dir) / 'report.html', html
                )
            _, browser_payload, embedded_data = self._browser_report_data(
                generator.output_dir
            )

            self.assertEqual(len(guide_data), 1005)
            self.assertEqual(guide_data[-1]['reference'], 'Single_1004')
            self.assertEqual(guide_data[-1]['target_sequence'], 'GUIDE1004NGG')
            self.assertEqual(len(embedded_data['guide_rna']), 1005)
            self.assertEqual(
                embedded_data['guide_rna'][-1]['reference'], 'Single_1004'
            )
            self.assertNotIn('guide_rna', browser_payload['data'])
            self.assertEqual(browser_payload['guides']['row_count'], 1005)

            javascript = (
                Path(generator.output_dir) / 'static' / 'js' / 'main.js'
            ).read_text(encoding='utf-8')
            self.assertNotIn('const rows = guideData.map', javascript)
            self.assertIn('this.setupGuideTableFilters(guideData);', javascript)

    def test_html_publish_failure_preserves_previous_payload_reference(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            final_guides_path = temp_dir / 'final_guides.tsv'
            self._write_final_guides(final_guides_path)
            results = self._results(
                config, final_annotation_guides=final_guides_path
            )

            first_generator = ReportGenerator(config)
            report_path = Path(first_generator.generate_html_report(results))
            first_html = report_path.read_bytes()
            first_core = Path(results['browser_payload']['core_path'])
            self.assertTrue(first_core.is_file())

            second_generator = ReportGenerator(config)
            original_atomic_write = second_generator._atomic_write_text

            def fail_report_commit(path, content):
                if Path(path).name == 'report.html':
                    raise OSError('simulated report commit failure')
                return original_atomic_write(path, content)

            with mock.patch.object(
                second_generator,
                '_atomic_write_text',
                side_effect=fail_report_commit,
            ):
                with self.assertRaisesRegex(
                    OSError, 'simulated report commit failure'
                ):
                    second_generator.generate_html_report(results)

            self.assertEqual(report_path.read_bytes(), first_html)
            self.assertTrue(first_core.is_file())
            failed_core = Path(results['browser_payload']['core_path'])
            self.assertTrue(failed_core.is_file())
            failed_marker = (
                Path(second_generator.web_data_dir)
                / 'report_payloads'
                / '.committed'
                / failed_core.parent.name
            )
            self.assertFalse(failed_marker.exists())

            final_generator = ReportGenerator(config)
            final_generator.generate_html_report(results)
            normal_generations = [
                path for path in (
                    Path(final_generator.web_data_dir) / 'report_payloads'
                ).iterdir()
                if path.is_dir() and path.name.startswith('normal-')
            ]
            self.assertEqual(len(normal_generations), 2)
            self.assertTrue(all(
                (generation / 'manifest.json').is_file()
                for generation in normal_generations
            ))
            self.assertFalse(failed_core.exists())
            self.assertTrue(first_core.is_file())
            _, payload, data = self._browser_report_data(
                final_generator.output_dir
            )
            self.assertEqual(payload['guides']['row_count'], 1)
            self.assertEqual(len(data['guide_rna']), 1)

    def test_multi_gene_contig_rows_are_counted_once_in_report_data(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            annotation_path = temp_dir / 'annotation_classified.txt'
            final_guides_path = temp_dir / 'annotation_with_guides_final.tsv'
            header = [
                'query_id', 'subject_id', 'identity', 'alignment_length',
                'evalue', 'bit_score', 'description', 'gene_name',
                'classification', 'unique_id',
            ]
            annotation_rows = [
                [
                    'contig_a', 'subject_multi', '99.0', '100', '1e-20',
                    '200', 'multi description', 'MULTI', 'Multi', 'Multi_1',
                ],
                [
                    'contig_b', 'subject_multi', '98.0', '90', '1e-18',
                    '180', 'multi description', 'MULTI', 'Multi', 'Multi_1',
                ],
                [
                    'contig_c', 'subject_single', '97.0', '80', '1e-16',
                    '160', 'single description', 'SINGLE', 'Single',
                    'Single_1',
                ],
            ]
            annotation_path.write_text(
                '\n'.join(
                    '\t'.join(row) for row in [header, *annotation_rows]
                ) + '\n',
                encoding='utf-8',
            )

            multi_guide_a = (
                'Multi_1/10/33/MULTI_GUIDE_A/1/0/1.0/+'
            )
            multi_guide_b = (
                'Multi_1/20/43/MULTI_GUIDE_B/1/0/1.0/-'
            )
            single_guide = (
                'Single_1/30/53/SINGLE_GUIDE/1/0/1.0/+'
            )
            final_rows = [
                [*annotation_rows[0], ','.join([
                    multi_guide_a, multi_guide_b, multi_guide_a,
                ])],
                [*annotation_rows[1], ','.join([
                    multi_guide_b, multi_guide_a,
                ])],
                [*annotation_rows[2], single_guide],
            ]
            final_guides_path.write_text(
                '\n'.join(
                    '\t'.join(row)
                    for row in [[*header, 'guide_targets'], *final_rows]
                ) + '\n',
                encoding='utf-8',
            )

            results = self._results(
                config,
                annotation_table=annotation_path,
                final_annotation_guides=final_guides_path,
            )
            generator = ReportGenerator(config)

            stats = generator._collect_statistics(results)
            gene_data = generator._load_gene_annotation_data(results)
            js_gene_data = generator._load_gene_annotation_data_for_js(
                results
            )
            guide_data = generator._load_guide_rna_data_for_js(results)

            self.assertEqual(stats['annotation'], {
                'multi_genes': 1,
                'single_genes': 1,
                'total_genes': 2,
            })
            self.assertEqual(stats['guide_design']['total_guides'], 3)
            self.assertEqual(
                stats['guide_design']['site_specific_guides'], 3
            )
            self.assertEqual(
                stats['guide_design']['site_specific_multi_target'], 2
            )
            self.assertEqual(
                stats['guide_design']['site_specific_single_target'], 1
            )

            self.assertEqual(len(gene_data), 2)
            self.assertEqual(len(js_gene_data), 2)
            multi_gene = next(
                gene for gene in gene_data
                if gene['unique_id'] == 'Multi_1'
            )
            self.assertEqual(multi_gene['query_id'], 'contig_a; contig_b')
            self.assertEqual(multi_gene['annotation_row_count'], 2)
            self.assertEqual(multi_gene['target_count'], 2)

            self.assertEqual(len(guide_data), 3)
            self.assertEqual(
                len({
                    (guide['unique_id'], guide['reference'], guide['start'],
                     guide['end'], guide['target_sequence'], guide['strand'])
                    for guide in guide_data
                }),
                3,
            )
            self.assertEqual(
                {
                    guide['query_id']
                    for guide in guide_data
                    if guide['unique_id'] == 'Multi_1'
                },
                {'contig_a; contig_b'},
            )

    def test_conflicting_rows_for_one_gene_group_fail_report_generation(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            config = self._config(temp_dir)
            annotation_path = temp_dir / 'annotation_classified.txt'
            annotation_path.write_text(
                '\t'.join([
                    'query_id', 'subject_id', 'identity',
                    'alignment_length', 'evalue', 'bit_score',
                    'description', 'gene_name', 'classification',
                    'unique_id',
                ])
                + '\n'
                + '\t'.join([
                    'contig_a', 'subject_a', '99.0', '100', '1e-20',
                    '200', 'description', 'GENE', 'Multi', 'Multi_1',
                ])
                + '\n'
                + '\t'.join([
                    'contig_b', 'subject_b', '98.0', '90', '1e-18',
                    '180', 'description', 'GENE', 'Multi', 'Multi_1',
                ])
                + '\n',
                encoding='utf-8',
            )
            generator = ReportGenerator(config)

            with self.assertRaisesRegex(
                ValueError, 'conflicting subject_id values'
            ):
                generator._collect_statistics(self._results(
                    config, annotation_table=annotation_path
                ))

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
