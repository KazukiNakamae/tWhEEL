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


class ReportDatabaseIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        output_dir = self.root / 'output'
        intermediate_dir = output_dir / 'intermediate'
        intermediate_dir.mkdir(parents=True)
        self.config = {
            'output': {
                'base_dir': str(output_dir),
                'intermediate': str(intermediate_dir),
            },
            'database': {'force_rebuild': False},
            'input': {},
            'crispr': {},
            'busco': {},
        }
        self.go_path = self.root / 'go_kfi_table.csv'
        self.generator = ReportGenerator(self.config)
        self.results = {
            'timestamp': '2026-07-16T00:00:00',
            'statistics': {},
            'config': self.config,
            'data': {'go_kfi_table': str(self.go_path)},
        }

    def tearDown(self):
        self.temporary_directory.cleanup()

    def _write_go_ids(self, *go_ids):
        rows = [
            'GO_ID,total_genes,genes_with_guides,total_guides,'
            'kfi_score,targetable_type'
        ]
        rows.extend(f'{go_id},1,1,1,1.0,Fully' for go_id in go_ids)
        self.go_path.write_text('\n'.join(rows) + '\n', encoding='utf-8')

    @staticmethod
    def _terms(prefix, *go_ids):
        return {
            go_id: {
                'name': f'{prefix} name {go_id}',
                'definition': f'{prefix} definition {go_id}',
            }
            for go_id in go_ids
        }

    def _prepare(self, api_terms):
        with mock.patch.object(
            self.generator,
            '_fetch_go_terms_from_api',
            return_value=api_terms,
        ):
            return self.generator.prepare_web_data(self.results)

    @staticmethod
    def _go_rows(db_path):
        with sqlite3.connect(db_path) as connection:
            return connection.execute(
                'SELECT go_id, go_name, go_definition '
                'FROM go_terms ORDER BY go_id'
            ).fetchall()

    def test_missing_go_row_forces_cache_rebuild(self):
        self._write_go_ids('GO:0000001')
        expected_terms = self._terms('cached', 'GO:0000001')
        db_path, _ = self._prepare(expected_terms)
        with sqlite3.connect(db_path) as connection:
            connection.execute(
                "DELETE FROM go_terms WHERE go_id = 'GO:0000001'"
            )

        with mock.patch.object(
            self.generator,
            '_fetch_go_terms_from_api',
            return_value=expected_terms,
        ), mock.patch.object(
            self.generator,
            '_rebuild_database',
            wraps=self.generator._rebuild_database,
        ) as rebuild_mock:
            self.generator.prepare_web_data(self.results)

        rebuild_mock.assert_called_once()
        self.assertEqual(
            self._go_rows(db_path),
            [(
                'GO:0000001',
                'cached name GO:0000001',
                'cached definition GO:0000001',
            )],
        )

    def test_extra_go_row_forces_cache_rebuild(self):
        self._write_go_ids('GO:0000001')
        expected_terms = self._terms('cached', 'GO:0000001')
        db_path, _ = self._prepare(expected_terms)
        with sqlite3.connect(db_path) as connection:
            connection.execute(
                "INSERT INTO go_terms VALUES ('GO:EXTRA', 'extra', 'extra')"
            )

        with mock.patch.object(
            self.generator,
            '_rebuild_database',
            wraps=self.generator._rebuild_database,
        ) as rebuild_mock:
            self.generator.prepare_web_data(self.results)

        rebuild_mock.assert_called_once()
        self.assertEqual(
            [row[0] for row in self._go_rows(db_path)],
            ['GO:0000001'],
        )

    def test_extra_metadata_key_forces_cache_rebuild(self):
        self._write_go_ids('GO:0000001')
        expected_terms = self._terms('cached', 'GO:0000001')
        db_path, _ = self._prepare(expected_terms)
        with sqlite3.connect(db_path) as connection:
            connection.execute(
                "INSERT INTO report_metadata VALUES ('unused', 'value')"
            )

        with mock.patch.object(
            self.generator,
            '_rebuild_database',
            wraps=self.generator._rebuild_database,
        ) as rebuild_mock:
            self.generator.prepare_web_data(self.results)

        rebuild_mock.assert_called_once()
        with sqlite3.connect(db_path) as connection:
            metadata_keys = {
                row[0] for row in connection.execute(
                    'SELECT key FROM report_metadata'
                )
            }
        self.assertEqual(
            metadata_keys,
            {'schema_version', 'input_fingerprint'},
        )

    def test_blank_go_row_forces_rebuild_and_recovers_content(self):
        self._write_go_ids('GO:0000001')
        db_path, _ = self._prepare({})
        self.assertEqual(
            self._go_rows(db_path),
            [('GO:0000001', '', '')],
        )
        recovered_terms = self._terms('recovered', 'GO:0000001')

        with mock.patch.object(
            self.generator,
            '_fetch_go_terms_from_api',
            return_value=recovered_terms,
        ), mock.patch.object(
            self.generator,
            '_rebuild_database',
            wraps=self.generator._rebuild_database,
        ) as rebuild_mock:
            self.generator.prepare_web_data(self.results)

        rebuild_mock.assert_called_once()
        self.assertEqual(
            self._go_rows(db_path),
            [(
                'GO:0000001',
                'recovered name GO:0000001',
                'recovered definition GO:0000001',
            )],
        )

    def test_force_refresh_keeps_cached_values_when_api_fails(self):
        self._write_go_ids('GO:0000001', 'GO:0000002')
        cached_terms = self._terms(
            'cached',
            'GO:0000001',
            'GO:0000002',
        )
        db_path, _ = self._prepare(cached_terms)
        original_rows = self._go_rows(db_path)
        self.generator.config['database']['force_rebuild'] = True

        with mock.patch.object(
            self.generator,
            '_fetch_go_terms_from_api',
            return_value={},
        ) as fetch_mock:
            self.generator.prepare_web_data(self.results)

        fetch_mock.assert_called_once_with([
            'GO:0000001',
            'GO:0000002',
        ])
        self.assertEqual(self._go_rows(db_path), original_rows)

    def test_lookup_writes_recovered_content_without_inserting_extra_ids(self):
        self._write_go_ids('GO:0000001')
        db_path, _ = self._prepare(
            self._terms('cached', 'GO:0000001')
        )
        with sqlite3.connect(db_path) as connection:
            connection.execute(
                "UPDATE go_terms SET go_name = '', go_definition = '' "
                "WHERE go_id = 'GO:0000001'"
            )
        api_terms = self._terms(
            'recovered',
            'GO:0000001',
            'GO:EXTRA',
        )

        with mock.patch.object(
            self.generator,
            '_fetch_go_terms_from_api',
            return_value=api_terms,
        ):
            recovered = self.generator._fetch_go_terms_from_database(
                ['GO:0000001', 'GO:EXTRA']
            )

        self.assertEqual(recovered, api_terms)
        self.assertEqual(
            self._go_rows(db_path),
            [(
                'GO:0000001',
                'recovered name GO:0000001',
                'recovered definition GO:0000001',
            )],
        )
        with mock.patch.object(
            self.generator,
            '_fetch_go_terms_from_api',
        ) as fetch_mock:
            self.generator._fetch_go_terms_from_database(['GO:0000001'])
        fetch_mock.assert_not_called()


if __name__ == '__main__':
    unittest.main()
