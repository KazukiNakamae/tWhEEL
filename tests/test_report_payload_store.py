import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from src.modules.report_payload_store import ReportPayloadStore


class ReportPayloadStoreTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temporary_directory.name)
        self.web_data_dir = self.temp_dir / 'web_data'

    def tearDown(self):
        self.temporary_directory.cleanup()

    def make_store(self, chunk_rows=5):
        return ReportPayloadStore(
            self.web_data_dir,
            guide_chunk_rows=chunk_rows,
        )

    @staticmethod
    def guide_rows(count, *, dangerous_value=None):
        for index in range(count):
            row = {
                'guide_id': f'guide-{index:04d}',
                'unique_id': f'Single_{index:010d}',
                'start': index,
            }
            if dangerous_value is not None:
                row['dangerous'] = dangerous_value
            yield row

    @staticmethod
    def extract_assignment(source, global_name, predicate=None):
        marker = f'window.{global_name}'
        decoder = json.JSONDecoder()
        offset = 0
        decoded_values = []
        while True:
            marker_index = source.find(marker, offset)
            if marker_index < 0:
                break
            equals_index = source.find('=', marker_index + len(marker))
            if equals_index < 0:
                break
            candidate = source[equals_index + 1:].lstrip()
            try:
                value, _ = decoder.raw_decode(candidate)
            except json.JSONDecodeError:
                offset = marker_index + len(marker)
                continue
            decoded_values.append(value)
            offset = equals_index + 1

        if predicate is not None:
            decoded_values = [value for value in decoded_values if predicate(value)]
        if len(decoded_values) != 1:
            raise AssertionError(
                f'Expected one JSON assignment to {global_name}, '
                f'found {len(decoded_values)}'
            )
        return decoded_values[0]

    @staticmethod
    def extract_variable_assignment(source, variable_name):
        marker = f'var {variable_name} = '
        marker_index = source.index(marker)
        candidate = source[marker_index + len(marker):].lstrip()
        value, _ = json.JSONDecoder().raw_decode(candidate)
        return value

    @classmethod
    def load_core_payload(cls, core_path):
        source = Path(core_path).read_text(encoding='utf-8')
        payload = cls.extract_variable_assignment(source, 'payload')
        if not isinstance(payload, dict) or 'generation_id' not in payload:
            raise AssertionError('Core JavaScript does not contain its payload')
        if 'window.__TWHEEL_REPORT_PAYLOAD__ = payload;' not in source:
            raise AssertionError('Core JavaScript does not publish its singleton')
        return source, payload

    @classmethod
    def load_chunk_payload(cls, chunk_path):
        source = Path(chunk_path).read_text(encoding='utf-8')
        payload = cls.extract_assignment(
            source,
            '__TWHEEL_REPORT_GUIDE_CHUNKS__',
            predicate=lambda value: isinstance(value, list),
        )
        return source, payload

    @staticmethod
    def manifest_file_path(manifest_path, relative_path):
        relative_path = Path(relative_path)
        if relative_path.is_absolute():
            return relative_path
        return Path(manifest_path).parent / relative_path

    def assert_manifest_file_integrity(
        self, manifest_path, relative_path, entry
    ):
        artifact_path = self.manifest_file_path(manifest_path, relative_path)
        self.assertTrue(artifact_path.is_file())
        content = artifact_path.read_bytes()
        self.assertEqual(entry['size'], len(content))
        self.assertEqual(
            entry['sha256'],
            hashlib.sha256(content).hexdigest(),
        )
        return artifact_path

    def test_publish_chunks_guides_and_records_complete_manifest(self):
        store = self.make_store(chunk_rows=5)
        core_data = {
            'timestamp': '2026-07-17T00:00:00',
            'statistics': {'guide_count': 12},
        }

        published = store.publish(
            'normal',
            core_data,
            self.guide_rows(12),
        )

        self.assertEqual(published['guide_count'], 12)
        self.assertEqual(len(published['chunks']), 3)
        self.assertRegex(published['generation_id'], r'^[0-9a-f]{32}$')
        self.assertTrue(Path(published['core_path']).is_file())
        self.assertTrue(Path(published['manifest_path']).is_file())
        self.assertEqual(
            published['core_relative_path'],
            'web_data/report_payloads/'
            f"normal-{published['generation_id']}/core.js",
        )

        core_source, core_payload = self.load_core_payload(
            published['core_path']
        )
        self.assertEqual(core_payload['namespace'], 'normal')
        self.assertEqual(
            core_payload['generation_id'], published['generation_id']
        )
        self.assertEqual(core_payload['data'], core_data)
        self.assertEqual(core_payload['guides']['row_count'], 12)
        self.assertEqual(core_payload['guides']['chunk_rows'], 5)
        self.assertEqual(
            core_payload['guides']['chunks'], published['chunks']
        )
        self.assertNotIn('guide-0000', core_source)
        self.assertNotIn('guide_rna', core_payload['data'])

        manifest_path = Path(published['manifest_path'])
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        self.assertEqual(manifest['namespace'], 'normal')
        self.assertEqual(
            manifest['generation_id'], published['generation_id']
        )
        self.assertEqual(manifest['guide_count'], 12)
        self.assertEqual(manifest['guide_chunk_rows'], 5)
        self.assertIn('core.js', manifest['files'])
        core_artifact = self.assert_manifest_file_integrity(
            manifest_path, 'core.js', manifest['files']['core.js']
        )
        self.assertEqual(core_artifact.resolve(), Path(
            published['core_path']
        ).resolve())

        chunk_entries = sorted(
            (
                (key, entry) for key, entry in manifest['files'].items()
                if key != 'core.js'
            ),
            key=lambda item: item[1]['index'],
        )
        self.assertEqual(
            [entry['row_count'] for _, entry in chunk_entries],
            [5, 5, 2],
        )
        self.assertEqual(
            [entry['registry_key'] for _, entry in chunk_entries],
            [chunk['registry_key'] for chunk in published['chunks']],
        )

        reconstructed_rows = []
        for relative_path, entry in chunk_entries:
            chunk_path = self.assert_manifest_file_integrity(
                manifest_path, relative_path, entry
            )
            _, rows = self.load_chunk_payload(chunk_path)
            self.assertEqual(len(rows), entry['row_count'])
            reconstructed_rows.extend(rows)
        self.assertEqual(len(reconstructed_rows), 12)
        self.assertEqual(reconstructed_rows[0]['guide_id'], 'guide-0000')
        self.assertEqual(reconstructed_rows[-1]['guide_id'], 'guide-0011')

    def test_javascript_assignments_preserve_dangerous_json_values(self):
        store = self.make_store(chunk_rows=2)
        dangerous = '</script>&\u2028\u2029'
        core_data = {'dangerous': dangerous}
        extra_globals = {
            '__TWHEEL_TEST_EXTRA__': {'dangerous': dangerous},
        }

        published = store.publish(
            'normal',
            core_data,
            self.guide_rows(1, dangerous_value=dangerous),
            extra_globals=extra_globals,
        )

        core_source, core_payload = self.load_core_payload(
            published['core_path']
        )
        self.assertEqual(core_payload['data']['dangerous'], dangerous)
        self.assertEqual(
            self.extract_variable_assignment(core_source, 'extraGlobals'),
            {'__TWHEEL_TEST_EXTRA__': {'dangerous': dangerous}},
        )

        manifest_path = Path(published['manifest_path'])
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        chunk_relative_path, chunk_entry = next(
            (key, entry) for key, entry in manifest['files'].items()
            if key != 'core.js'
        )
        chunk_path = self.manifest_file_path(
            manifest_path, chunk_relative_path
        )
        chunk_source, chunk_rows = self.load_chunk_payload(chunk_path)
        self.assertEqual(chunk_rows[0]['dangerous'], dangerous)

        for source in (core_source, chunk_source):
            self.assertNotIn('</script>', source.lower())
            self.assertNotIn(dangerous, source)
            self.assertIn('\\u003c/script\\u003e\\u0026', source)
            self.assertNotIn('\u2028', source.replace('\\u2028', ''))
            self.assertNotIn('\u2029', source.replace('\\u2029', ''))

    def test_empty_guides_publish_core_without_chunks(self):
        store = self.make_store(chunk_rows=5)

        published = store.publish(
            'normal',
            {'statistics': {'guide_count': 0}},
            iter(()),
        )

        self.assertEqual(published['guide_count'], 0)
        self.assertEqual(published['chunks'], [])
        _, core_payload = self.load_core_payload(published['core_path'])
        self.assertEqual(core_payload['guides']['row_count'], 0)
        self.assertEqual(core_payload['guides']['chunks'], [])

        manifest = json.loads(
            Path(published['manifest_path']).read_text(encoding='utf-8')
        )
        self.assertEqual(manifest['guide_count'], 0)
        self.assertEqual(set(manifest['files']), {'core.js'})

    def test_guide_iterator_failure_leaves_no_final_or_temporary_generation(self):
        store = self.make_store(chunk_rows=2)

        def failing_rows():
            yield {'guide_id': 'first'}
            yield {'guide_id': 'second'}
            raise RuntimeError('simulated guide iterator failure')

        with self.assertRaisesRegex(
            RuntimeError, 'simulated guide iterator failure'
        ):
            store.publish('normal', {'status': 'new'}, failing_rows())

        payload_root = self.web_data_dir / 'report_payloads'
        if payload_root.exists():
            self.assertEqual(
                [
                    path for path in payload_root.iterdir()
                    if path.name.startswith('normal-')
                    or path.name.startswith('.normal-')
                ],
                [],
            )

    def test_nonfinite_numbers_fail_without_publishing_javascript(self):
        store = self.make_store(chunk_rows=2)

        with self.assertRaisesRegex(ValueError, 'Out of range float values'):
            store.publish(
                'normal',
                {'invalid_score': float('nan')},
                self.guide_rows(1),
            )

        payload_root = self.web_data_dir / 'report_payloads'
        self.assertEqual(
            [
                path for path in payload_root.iterdir()
                if path.is_dir()
            ],
            [],
        )

    def test_cleanup_normal_namespace_preserves_pathway_generations(self):
        store = self.make_store(chunk_rows=5)
        first_normal = store.publish(
            'normal', {'version': 1}, self.guide_rows(1)
        )
        pathway = store.publish(
            'pathway', {'version': 1}, self.guide_rows(1)
        )
        second_normal = store.publish(
            'normal', {'version': 2}, self.guide_rows(1)
        )

        first_normal_dir = Path(first_normal['core_path']).parent
        pathway_dir = Path(pathway['core_path']).parent
        second_normal_dir = Path(second_normal['core_path']).parent
        self.assertTrue(first_normal_dir.is_dir())
        self.assertTrue(pathway_dir.is_dir())
        self.assertTrue(second_normal_dir.is_dir())

        with store.publication_lock():
            store.cleanup_except('normal', second_normal['generation_id'])

        self.assertFalse(first_normal_dir.exists())
        self.assertTrue(second_normal_dir.is_dir())
        self.assertTrue(pathway_dir.is_dir())
        self.assertTrue(Path(pathway['manifest_path']).is_file())

    def test_orphan_cleanup_preserves_committed_immutable_generation(self):
        store = self.make_store(chunk_rows=5)
        with store.publication_lock():
            committed = store.publish(
                'normal', {'version': 1}, self.guide_rows(1)
            )
            committed_dir = Path(committed['core_path']).parent
            files_before_commit = sorted(
                path.relative_to(committed_dir).as_posix()
                for path in committed_dir.rglob('*')
            )
            store.mark_committed('normal', committed['generation_id'])
            orphan = store.publish(
                'normal', {'version': 2}, self.guide_rows(1)
            )
            current = store.publish(
                'normal', {'version': 3}, self.guide_rows(1)
            )
            store.cleanup_uncommitted_locked(
                'normal', current['generation_id']
            )

        self.assertTrue(committed_dir.is_dir())
        self.assertFalse(Path(orphan['core_path']).parent.exists())
        self.assertTrue(Path(current['core_path']).parent.is_dir())
        self.assertEqual(
            files_before_commit,
            sorted(
                path.relative_to(committed_dir).as_posix()
                for path in committed_dir.rglob('*')
            ),
        )
        marker = (
            store.committed_root
            / f"normal-{committed['generation_id']}"
        )
        self.assertTrue(marker.is_file())


if __name__ == '__main__':
    unittest.main()
