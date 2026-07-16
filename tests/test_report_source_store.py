import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from src.modules.report_source_store import ReportSourceStore


class ReportSourceStoreTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temporary_directory.name)
        self.web_data_dir = self.temp_dir / 'web_data'
        self.store = ReportSourceStore(self.web_data_dir)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def _source(self, name, content):
        path = self.temp_dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding='utf-8')
        return path

    def _manifest(self):
        return json.loads(self.store.manifest_path.read_text(encoding='utf-8'))

    def test_publish_creates_valid_version_two_snapshot(self):
        table = self._source('inputs/table.tsv', 'id\tvalue\n1\talpha\n')
        image = self._source('inputs/image.bin', 'image payload')

        published = self.store.publish_namespace('normal', {
            'table': (table, 'tables/table.tsv'),
            'image': (image, 'images/image.bin'),
        })

        self.assertEqual(set(published), {'table', 'image'})
        self.assertEqual(self.store.load_namespace('normal'), published)
        table_snapshot = Path(published['table'])
        image_snapshot = Path(published['image'])
        self.assertFalse(os.path.samefile(table, table_snapshot))
        self.assertFalse(os.path.samefile(image, image_snapshot))
        self.assertEqual(
            table_snapshot.read_text(encoding='utf-8'),
            table.read_text(encoding='utf-8'),
        )

        manifest = self._manifest()
        self.assertEqual(manifest['version'], 2)
        self.assertEqual(
            set(manifest['sources']),
            {'normal.table', 'normal.image'},
        )
        table_entry = manifest['sources']['normal.table']
        self.assertFalse(Path(table_entry['path']).is_absolute())
        self.assertEqual(table_entry['size'], table.stat().st_size)
        self.assertEqual(
            table_entry['sha256'],
            hashlib.sha256(table.read_bytes()).hexdigest(),
        )

    def test_overwriting_source_does_not_change_snapshot(self):
        source = self._source('large.tsv', 'original snapshot\n')
        published = self.store.publish_namespace('normal', {
            'large_table': (source, 'large.tsv'),
        })
        snapshot = Path(published['large_table'])

        source.write_text('replacement data\n', encoding='utf-8')

        self.assertEqual(snapshot.read_text(encoding='utf-8'), 'original snapshot\n')
        self.assertEqual(
            self.store.load_namespace('normal'),
            {'large_table': str(snapshot.resolve())},
        )

    @unittest.skipUnless(sys.platform == 'darwin', 'requires macOS clonefile cp')
    def test_macos_publish_uses_clone_copy_with_an_independent_inode(self):
        source = self._source('clone-source.tsv', 'clone payload\n')
        original_run = subprocess.run

        with mock.patch(
            'src.modules.report_source_store.subprocess.run',
            wraps=original_run,
        ) as run_mock:
            published = self.store.publish_namespace('normal', {
                'table': (source, 'table.tsv'),
            })

        command = run_mock.call_args.args[0]
        self.assertEqual(command[:3], ['cp', '-c', '-p'])
        self.assertFalse(os.path.samefile(source, published['table']))

    def test_replace_prunes_namespace_and_preserves_other_namespace(self):
        normal_a = self._source('normal-a.tsv', 'normal a v1\n')
        normal_b = self._source('normal-b.tsv', 'normal b v1\n')
        pathway = self._source('pathway.tsv', 'pathway v1\n')
        first_normal = self.store.publish_namespace('normal', {
            'a': (normal_a, 'a.tsv'),
            'b': (normal_b, 'b.tsv'),
        })
        published_pathway = self.store.publish_namespace('pathway', {
            'mapping': (pathway, 'mapping.tsv'),
        })
        old_normal_snapshot = Path(first_normal['a']).parent
        pathway_snapshot = Path(published_pathway['mapping'])

        normal_a_v2 = self._source('normal-a-v2.tsv', 'normal a v2\n')
        second_normal = self.store.publish_namespace('normal', {
            'a': (normal_a_v2, 'a.tsv'),
        })

        self.assertEqual(set(second_normal), {'a'})
        self.assertFalse(old_normal_snapshot.exists())
        self.assertTrue(pathway_snapshot.is_file())
        self.assertEqual(
            self.store.load_namespace('pathway'),
            published_pathway,
        )
        self.assertEqual(
            set(self._manifest()['sources']),
            {'normal.a', 'pathway.mapping'},
        )

        second_normal_snapshot = Path(second_normal['a']).parent
        self.assertEqual(self.store.publish_namespace('normal', {}), {})
        self.assertEqual(self.store.load_namespace('normal'), {})
        self.assertFalse(second_normal_snapshot.exists())
        self.assertTrue(pathway_snapshot.is_file())
        self.assertEqual(
            set(self._manifest()['sources']),
            {'pathway.mapping'},
        )

    def test_manifest_publish_failure_preserves_old_manifest_and_snapshot(self):
        old_source = self._source('old.tsv', 'old generation\n')
        old_published = self.store.publish_namespace('normal', {
            'table': (old_source, 'table.tsv'),
        })
        old_snapshot = Path(old_published['table'])
        old_manifest = self.store.manifest_path.read_bytes()
        new_source = self._source('new.tsv', 'new generation\n')

        with mock.patch.object(
            self.store,
            '_write_manifest_atomic',
            side_effect=OSError('manifest publish failed'),
        ):
            with self.assertRaisesRegex(OSError, 'manifest publish failed'):
                self.store.publish_namespace('normal', {
                    'table': (new_source, 'table.tsv'),
                })

        self.assertEqual(self.store.manifest_path.read_bytes(), old_manifest)
        self.assertTrue(old_snapshot.is_file())
        self.assertEqual(old_snapshot.read_text(encoding='utf-8'), 'old generation\n')
        self.assertEqual(self.store.load_namespace('normal'), old_published)
        snapshot_directories = [
            path for path in self.store.snapshots_root.iterdir()
            if path.is_dir() and not path.name.startswith('.')
        ]
        self.assertEqual(
            [path.resolve() for path in snapshot_directories],
            [old_snapshot.parent.resolve()],
        )

    def test_copy_failure_preserves_old_manifest_and_snapshot(self):
        old_source = self._source('old.tsv', 'old generation\n')
        old_published = self.store.publish_namespace('normal', {
            'table': (old_source, 'table.tsv'),
        })
        old_snapshot = Path(old_published['table'])
        old_manifest = self.store.manifest_path.read_bytes()
        first_source = self._source('first.tsv', 'first new file\n')
        second_source = self._source('second.tsv', 'second new file\n')
        original_copy = self.store._copy_independent

        def fail_on_second_copy(source, destination):
            if destination.name == 'second.tsv':
                raise OSError('snapshot copy failed')
            original_copy(source, destination)

        with mock.patch.object(
            self.store,
            '_copy_independent',
            side_effect=fail_on_second_copy,
        ):
            with self.assertRaisesRegex(OSError, 'snapshot copy failed'):
                self.store.publish_namespace('normal', {
                    'first': (first_source, 'first.tsv'),
                    'second': (second_source, 'second.tsv'),
                })

        self.assertEqual(self.store.manifest_path.read_bytes(), old_manifest)
        self.assertTrue(old_snapshot.is_file())
        self.assertEqual(self.store.load_namespace('normal'), old_published)
        temporary_snapshots = [
            path for path in self.store.snapshots_root.iterdir()
            if path.name.startswith('.')
        ]
        self.assertEqual(temporary_snapshots, [])

    def test_load_detects_same_size_content_tampering(self):
        source = self._source('source.tsv', 'original\n')
        published = self.store.publish_namespace('normal', {
            'table': (source, 'table.tsv'),
        })
        snapshot = Path(published['table'])

        snapshot.write_text('tampered\n', encoding='utf-8')

        self.assertEqual(snapshot.stat().st_size, source.stat().st_size)
        with self.assertRaisesRegex(RuntimeError, 'checksum mismatch'):
            self.store.load_namespace('normal')

    def test_load_rejects_manifest_path_outside_snapshot_storage(self):
        source = self._source('source.tsv', 'payload\n')
        self.store.publish_namespace('normal', {
            'table': (source, 'table.tsv'),
        })
        manifest = self._manifest()
        manifest['sources']['normal.table']['path'] = '../source.tsv'
        self.store.manifest_path.write_text(
            json.dumps(manifest),
            encoding='utf-8',
        )

        with self.assertRaisesRegex(RuntimeError, 'escapes web_data'):
            self.store.load_namespace('normal')


if __name__ == '__main__':
    unittest.main()
