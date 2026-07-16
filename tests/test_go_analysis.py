import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from src.modules.gene_guide_normalization import (
    publish_normalized_gene_guide_generation,
)
from src.modules.go_analysis import GOAnalysis


class GoAnnotationTests(unittest.TestCase):
    def test_logs_uniprot_network_requirement_before_annotation(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            analysis = GOAnalysis({
                'output': {
                    'intermediate': temp_dir_name,
                },
                'go_analysis': {
                    'add_go_image': 'add-go-columns:test',
                },
            })
            completed = subprocess.CompletedProcess(
                args=[],
                returncode=0,
                stdout='',
                stderr='',
            )

            with self.assertLogs(
                'src.modules.go_analysis',
                level='INFO',
            ) as captured_logs:
                with mock.patch(
                    'src.modules.go_analysis.subprocess.run',
                    return_value=completed,
                ):
                    analysis.add_go_annotations('input.tsv')

        self.assertIn(
            'GO annotation requires network access to the UniProt REST API: '
            'https://rest.uniprot.org',
            '\n'.join(captured_logs.output),
        )

    def test_resolves_published_normalized_genes_before_legacy_table(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            annotation_path = temp_dir / 'annotation_classified.txt'
            annotation_path.write_text(
                'query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\t'
                'description\tgene_name\tclassification\tunique_id\n'
                'query-1\tsp|P04637|P53_HUMAN\t99\t100\t1e-20\t200\t'
                'description\tTP53\tSingle\tSingle_1\n',
                encoding='utf-8',
            )
            single_bed = temp_dir / 'single.bed'
            single_bed.write_text(
                'Single_1\t1\t23\tACGT/1/0\t1\t+\n',
                encoding='utf-8',
            )
            published = publish_normalized_gene_guide_generation(
                annotation_path,
                temp_dir,
                single_bed_path=single_bed,
            )
            (temp_dir / 'annotation_with_guides_final.tsv').write_text(
                'stale legacy table\n', encoding='utf-8'
            )
            analysis = GOAnalysis({
                'output': {'intermediate': str(temp_dir)},
            })

            input_path, input_format = (
                analysis._resolve_gene_annotation_input()
            )

            self.assertEqual(input_format, 'normalized')
            self.assertEqual(
                Path(input_path), Path(published['normalized_genes'])
            )

    def test_invalid_normalized_pointer_does_not_fall_back_to_legacy(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            (temp_dir / 'annotation_with_guides_final.tsv').write_text(
                'legacy\n', encoding='utf-8'
            )
            (temp_dir / 'normalized_gene_guide_current.json').write_text(
                '{invalid json', encoding='utf-8'
            )
            analysis = GOAnalysis({
                'output': {'intermediate': str(temp_dir)},
            })

            with self.assertRaisesRegex(
                ValueError, 'Cannot read normalized generation pointer'
            ):
                analysis._resolve_gene_annotation_input()

    def test_broken_pointer_symlink_does_not_fall_back_to_legacy(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            (temp_dir / 'annotation_with_guides_final.tsv').write_text(
                'legacy\n', encoding='utf-8'
            )
            (temp_dir / 'normalized_gene_guide_current.json').symlink_to(
                'missing-pointer.json'
            )
            analysis = GOAnalysis({
                'output': {'intermediate': str(temp_dir)},
            })

            with self.assertRaisesRegex(
                ValueError, 'Invalid normalized generation pointer'
            ):
                analysis._resolve_gene_annotation_input()

    def test_legacy_input_remains_available_without_current_pointer(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            legacy_path = temp_dir / 'annotation_with_guides_final.tsv'
            legacy_path.write_text('legacy\n', encoding='utf-8')
            analysis = GOAnalysis({
                'output': {'intermediate': str(temp_dir)},
            })

            input_path, input_format = (
                analysis._resolve_gene_annotation_input()
            )

            self.assertEqual(input_format, 'legacy')
            self.assertEqual(Path(input_path), legacy_path)

    def test_missing_helper_image_build_does_not_require_sudo(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            analysis = GOAnalysis({
                'output': {'intermediate': temp_dir_name},
                'go_analysis': {
                    'add_go_image': 'add-go-columns:v1.2-test',
                    'compute_kfi_image': 'compute-go-kfi:existing',
                    'count_parent_image': 'count-parent:existing',
                },
            })
            completed = subprocess.CompletedProcess(
                args=[], returncode=0, stdout='', stderr=''
            )

            with mock.patch.object(
                analysis,
                '_check_docker_image_exists',
                side_effect=[False, True, True],
            ), mock.patch(
                'src.modules.go_analysis.subprocess.run',
                return_value=completed,
            ) as run_mock:
                analysis.build_docker_images()

            command = run_mock.call_args.args[0]
            self.assertEqual(command[:2], ['docker', 'build'])
            self.assertNotIn('sudo', command)

    def test_cached_legacy_helpers_are_rebuilt_for_normalized_input(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            analysis = GOAnalysis({
                'output': {'intermediate': temp_dir_name},
                'go_analysis': {
                    'add_go_image': 'custom-add-go:legacy',
                    'compute_kfi_image': 'custom-compute-kfi:legacy',
                    'count_parent_image': 'count-parent:existing',
                },
            })
            completed = subprocess.CompletedProcess(
                args=[], returncode=0, stdout='', stderr=''
            )

            with mock.patch.object(
                analysis,
                '_check_docker_image_exists',
                return_value=True,
            ), mock.patch.object(
                analysis,
                '_image_supports_normalized_gene_guides',
                side_effect=[False, True, False, True],
            ) as capability_mock, mock.patch(
                'src.modules.go_analysis.subprocess.run',
                return_value=completed,
            ) as run_mock:
                analysis.build_docker_images(
                    require_normalized_gene_guides=True
                )

            self.assertEqual(capability_mock.call_count, 4)
            self.assertEqual(run_mock.call_count, 2)
            for call in run_mock.call_args_list:
                self.assertEqual(call.args[0][:2], ['docker', 'build'])


class GoParentTermAnalysisTests(unittest.TestCase):
    @staticmethod
    def _config(intermediate_dir):
        return {
            'output': {
                'intermediate': str(intermediate_dir),
            },
            'go_analysis': {
                'count_parent_image': 'count-go-parent:test',
            },
        }

    def test_parent_term_docker_command_does_not_request_a_tty(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            fully_targetable_path = temp_dir / 'fully_targetable.csv'
            fully_targetable_path.write_text(
                'GO_ID,total_genes\nGO:0000001,1\n',
                encoding='utf-8',
            )
            analysis = GOAnalysis(self._config(temp_dir))

            completed = subprocess.CompletedProcess(
                args=[],
                returncode=0,
                stdout='completed',
                stderr='',
            )
            with mock.patch(
                'src.modules.go_analysis.subprocess.run',
                return_value=completed,
            ) as run_mock:
                csv_path, png_path = analysis.count_go_parent_terms(
                    str(fully_targetable_path)
                )

            command = run_mock.call_args.args[0]
            self.assertEqual(command[:3], ['docker', 'run', '--rm'])
            self.assertNotIn('-i', command)
            self.assertNotIn('-t', command)
            self.assertNotIn('-it', command)
            self.assertEqual(
                run_mock.call_args.kwargs,
                {'capture_output': True, 'text': True},
            )
            self.assertEqual(csv_path, str(temp_dir / 'go_parent_counts.csv'))
            self.assertEqual(png_path, str(temp_dir / 'go_parent_counts.png'))

    def test_parent_term_container_failure_is_propagated(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            fully_targetable_path = temp_dir / 'fully_targetable.csv'
            fully_targetable_path.write_text(
                'GO_ID,total_genes\nGO:0000001,1\n',
                encoding='utf-8',
            )
            analysis = GOAnalysis(self._config(temp_dir))
            failed = subprocess.CompletedProcess(
                args=[],
                returncode=1,
                stdout='',
                stderr='parent term analysis error',
            )

            with mock.patch(
                'src.modules.go_analysis.subprocess.run',
                return_value=failed,
            ):
                with self.assertRaisesRegex(
                    RuntimeError,
                    'GO parent term analysis failed',
                ):
                    analysis.count_go_parent_terms(
                        str(fully_targetable_path)
                    )

    def test_empty_input_skips_docker_in_non_tty_execution(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            fully_targetable_path = temp_dir / 'fully_targetable.csv'
            fully_targetable_path.write_text(
                'GO_ID,total_genes\n',
                encoding='utf-8',
            )
            analysis = GOAnalysis(self._config(temp_dir))

            with mock.patch(
                'src.modules.go_analysis.subprocess.run'
            ) as run_mock:
                csv_path, message_path = analysis.count_go_parent_terms(
                    str(fully_targetable_path)
                )

            run_mock.assert_not_called()
            self.assertEqual(
                Path(csv_path).read_text(encoding='utf-8'),
                'go_term,count,parent_term\n',
            )
            self.assertTrue(Path(message_path).is_file())


class GoReverseLookupTests(unittest.TestCase):
    @staticmethod
    def _config(intermediate_dir):
        return {
            'output': {
                'intermediate': str(intermediate_dir),
            },
        }

    def test_uses_go_all_once_and_writes_only_unique_gene_ids(self):
        header = [
            'query_id',
            'subject_id',
            'identity',
            'length',
            'evalue',
            'bit_score',
            'description',
            'gene_name',
            'classification',
            'unique_id',
            'guide_targets',
            'GO_BP',
            'GO_MF',
            'GO_CC',
            'GO_ALL',
        ]

        def make_row(**values):
            return '\t'.join(values.get(column, '') for column in header)

        rows = [
            make_row(
                query_id='query1',
                unique_id='Single_1',
                description='description that must not be copied',
                guide_targets='guide_a,guide_b',
                GO_BP='GO:0000001;GO:0000002',
                GO_ALL='GO:0000001;GO:0000002;GO:0000001',
            ),
            make_row(
                query_id='query1_duplicate',
                unique_id='Single_1',
                guide_targets='different_guide',
                GO_BP='GO:0000001',
                GO_ALL='GO:0000001',
            ),
            make_row(
                query_id='query2',
                unique_id='Single_2',
                guide_targets='guide_c',
                GO_BP='GO:0000001',
                GO_ALL='GO:0000001',
            ),
            make_row(
                query_id='query3',
                unique_id='Single_3',
                guide_targets='guide_d',
                GO_MF='GO:0000003',
                GO_ALL='',
            ),
        ]

        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'annotation_with_guides_go.tsv'
            input_path.write_text(
                '\n'.join(['\t'.join(header)] + rows) + '\n',
                encoding='utf-8',
            )
            analysis = GOAnalysis(self._config(temp_dir))

            output_path = Path(analysis.create_go_reverse_lookup(input_path))
            output_lines = output_path.read_text(encoding='utf-8').splitlines()

        self.assertEqual(
            output_lines,
            [
                'go_term\tgene_count\tgenes',
                'GO:0000001\t2\tSingle_1;Single_2',
                'GO:0000002\t1\tSingle_1',
                'GO:0000003\t1\tSingle_3',
            ],
        )
        output_text = '\n'.join(output_lines)
        self.assertNotIn('guide_a', output_text)
        self.assertNotIn('description that must not be copied', output_text)
        self.assertEqual(output_text.count('Single_1'), 2)


if __name__ == '__main__':
    unittest.main()
