import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from src.modules.annotation import Annotation
from src.modules.sgd_annotation import Sgd_annotation


class FastaSearchCacheTests(unittest.TestCase):
    @staticmethod
    def _config(intermediate_dir):
        return {
            'output': {
                'intermediate': str(intermediate_dir),
            },
            'fasta_search': {
                'docker_image': 'fasta:test',
                'threads': 2,
                'evalue': 1e-5,
            },
        }

    @staticmethod
    def _successful_search(contents):
        def run(cmd, stdout=None, stderr=None, text=None):
            stdout.write(contents)
            return subprocess.CompletedProcess(cmd, 0, stderr='')

        return run

    def test_annotation_modules_keep_distinct_validated_caches(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            query_path = temp_dir / 'query.fasta'
            uniprot_path = temp_dir / 'uniprot.fasta'
            sgd_path = temp_dir / 'sgd.fasta'
            legacy_path = temp_dir / 'ggsearch_results.txt'
            query_path.write_text('>query\nMPEPTIDE\n', encoding='utf-8')
            uniprot_path.write_text('>uniprot\nMPEPTIDE\n', encoding='utf-8')
            sgd_path.write_text('>sgd\nMPEPTIDE\n', encoding='utf-8')
            legacy_path.write_text('legacy shared result\n', encoding='utf-8')

            annotator = Annotation(self._config(temp_dir))
            sgd_annotator = Sgd_annotation(self._config(temp_dir))

            def database_specific_search(cmd, stdout=None, stderr=None, text=None):
                stdout.write(f'result for {cmd[-1]}\n')
                return subprocess.CompletedProcess(cmd, 0, stderr='')

            with mock.patch(
                'src.modules.fasta_search_cache.subprocess.run',
                side_effect=database_specific_search,
            ) as run_mock:
                annotation_result = annotator.run_fasta_search(
                    str(query_path),
                    str(uniprot_path),
                )
                sgd_result = sgd_annotator.run_fasta_search(
                    str(query_path),
                    str(sgd_path),
                )

                self.assertEqual(run_mock.call_count, 2)
                self.assertEqual(
                    Path(annotation_result).name,
                    'ggsearch_results_annotation.txt',
                )
                self.assertEqual(
                    Path(sgd_result).name,
                    'ggsearch_results_sgd.txt',
                )
                self.assertNotEqual(annotation_result, sgd_result)
                self.assertEqual(
                    Path(annotation_result).read_text(encoding='utf-8'),
                    'result for uniprot.fasta\n',
                )
                self.assertEqual(
                    Path(sgd_result).read_text(encoding='utf-8'),
                    'result for sgd.fasta\n',
                )

                self.assertEqual(
                    annotator.run_fasta_search(str(query_path), str(uniprot_path)),
                    annotation_result,
                )
                self.assertEqual(
                    sgd_annotator.run_fasta_search(str(query_path), str(sgd_path)),
                    sgd_result,
                )
                self.assertEqual(run_mock.call_count, 2)

            annotation_metadata = json.loads(
                (temp_dir / 'ggsearch_results_annotation.meta.json').read_text(
                    encoding='utf-8'
                )
            )
            sgd_metadata = json.loads(
                (temp_dir / 'ggsearch_results_sgd.meta.json').read_text(
                    encoding='utf-8'
                )
            )
            self.assertEqual(annotation_metadata['namespace'], 'annotation')
            self.assertEqual(sgd_metadata['namespace'], 'sgd')
            self.assertEqual(
                annotation_metadata['database']['path'],
                str(uniprot_path.resolve()),
            )
            self.assertEqual(
                sgd_metadata['database']['path'],
                str(sgd_path.resolve()),
            )
            self.assertEqual(
                legacy_path.read_text(encoding='utf-8'),
                'legacy shared result\n',
            )

    def test_database_change_invalidates_existing_module_cache(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            query_path = temp_dir / 'query.fasta'
            database_path = temp_dir / 'database.fasta'
            query_path.write_text('>query\nMPEPTIDE\n', encoding='utf-8')
            database_path.write_text('>db\nMPEPTIDE\n', encoding='utf-8')
            annotator = Annotation(self._config(temp_dir))

            search_results = iter(['first result\n', 'updated result\n'])

            def changing_search(cmd, stdout=None, stderr=None, text=None):
                stdout.write(next(search_results))
                return subprocess.CompletedProcess(cmd, 0, stderr='')

            with mock.patch(
                'src.modules.fasta_search_cache.subprocess.run',
                side_effect=changing_search,
            ) as run_mock:
                result_path = annotator.run_fasta_search(
                    str(query_path),
                    str(database_path),
                )
                database_path.write_text(
                    '>db replacement\nMPEPTIDEWITHCHANGE\n',
                    encoding='utf-8',
                )
                rerun_path = annotator.run_fasta_search(
                    str(query_path),
                    str(database_path),
                )

            self.assertEqual(run_mock.call_count, 2)
            self.assertEqual(rerun_path, result_path)
            self.assertEqual(
                Path(rerun_path).read_text(encoding='utf-8'),
                'updated result\n',
            )

    def test_failed_search_does_not_publish_partial_cache(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            query_path = temp_dir / 'query.fasta'
            database_path = temp_dir / 'database.fasta'
            query_path.write_text('>query\nMPEPTIDE\n', encoding='utf-8')
            database_path.write_text('>db\nMPEPTIDE\n', encoding='utf-8')
            annotator = Annotation(self._config(temp_dir))

            def failed_search(cmd, stdout=None, stderr=None, text=None):
                stdout.write('partial result\n')
                return subprocess.CompletedProcess(cmd, 7, stderr='ggsearch error')

            with mock.patch(
                'src.modules.fasta_search_cache.subprocess.run',
                side_effect=failed_search,
            ):
                with self.assertRaisesRegex(RuntimeError, 'exit code 7'):
                    annotator.run_fasta_search(
                        str(query_path),
                        str(database_path),
                    )

            result_path = temp_dir / 'ggsearch_results_annotation.txt'
            metadata_path = temp_dir / 'ggsearch_results_annotation.meta.json'
            self.assertFalse(result_path.exists())
            self.assertFalse(metadata_path.exists())
            self.assertEqual(list(temp_dir.glob('*.tmp')), [])

            with mock.patch(
                'src.modules.fasta_search_cache.subprocess.run',
                side_effect=self._successful_search('complete result\n'),
            ) as run_mock:
                published_path = annotator.run_fasta_search(
                    str(query_path),
                    str(database_path),
                )

            self.assertEqual(run_mock.call_count, 1)
            self.assertEqual(
                Path(published_path).read_text(encoding='utf-8'),
                'complete result\n',
            )
            self.assertTrue(metadata_path.is_file())


if __name__ == '__main__':
    unittest.main()
