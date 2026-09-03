import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from src.modules.docker_workspace import DockerWorkspace
from src.modules.gene_prediction import GenePrediction
from src.modules.go_analysis import GOAnalysis
from src.modules.guide_design import GuideDesign


class DockerWorkspaceTests(unittest.TestCase):
    def test_absolute_and_relative_host_paths_map_to_same_container_path(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            root = Path(temp_dir_name) / 'output with spaces' / 'intermediate'
            root.mkdir(parents=True)
            host_path = root / 'nested' / 'result.tsv'
            workspace = DockerWorkspace(root, '/work')
            relative_root_workspace = DockerWorkspace(
                os.path.relpath(root), '/work'
            )

            self.assertEqual(
                workspace.container_path(host_path),
                '/work/nested/result.tsv',
            )
            self.assertEqual(
                workspace.container_path(os.path.relpath(host_path)),
                '/work/nested/result.tsv',
            )
            self.assertEqual(
                relative_root_workspace.container_path(host_path),
                '/work/nested/result.tsv',
            )
            self.assertEqual(
                workspace.volume_spec,
                f'{root.resolve()}:/work',
            )

    def test_path_outside_workspace_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            root = Path(temp_dir_name) / 'intermediate'
            root.mkdir()
            workspace = DockerWorkspace(root, '/work')

            with self.assertRaisesRegex(ValueError, 'outside the mounted'):
                workspace.container_path(root.parent / 'outside.tsv')

    def test_symlink_escape_from_workspace_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            parent = Path(temp_dir_name)
            root = parent / 'intermediate'
            outside = parent / 'outside'
            root.mkdir()
            outside.mkdir()
            (root / 'escape').symlink_to(outside, target_is_directory=True)
            workspace = DockerWorkspace(root, '/work')

            with self.assertRaisesRegex(ValueError, 'outside the mounted'):
                workspace.container_path(root / 'escape' / 'result.tsv')


class AbsoluteOutputDockerCommandTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir_context = tempfile.TemporaryDirectory()
        self.output_dir = (
            Path(self.temp_dir_context.name) / 'absolute output root'
        )
        self.intermediate_dir = self.output_dir / 'intermediate'
        self.intermediate_dir.mkdir(parents=True)
        self.completed = subprocess.CompletedProcess(
            args=[], returncode=0, stdout='completed', stderr=''
        )

    def tearDown(self):
        self.temp_dir_context.cleanup()

    def go_config(self):
        return {
            'output': {'intermediate': str(self.intermediate_dir)},
            'go_analysis': {
                'add_go_image': 'add-go:test',
                'compute_kfi_image': 'compute-kfi:test',
                'count_parent_image': 'count-parent:test',
            },
        }

    def assert_workspace_command(self, command, container_root):
        mount_index = command.index('-v') + 1
        self.assertEqual(
            command[mount_index],
            f'{self.intermediate_dir.resolve()}:{container_root}',
        )
        for argument in command[mount_index + 1:]:
            self.assertNotIn(str(self.intermediate_dir.resolve()), argument)
            self.assertNotIn(f'{container_root}//', argument)

    def test_go_commands_use_intermediate_workspace_for_all_io(self):
        input_path = self.intermediate_dir / 'annotation.tsv'
        input_path.write_text('header\n', encoding='utf-8')
        fully_path = self.intermediate_dir / 'fully_targetable.csv'
        fully_path.write_text(
            'GO_ID,total_genes\nGO:0000001,1\n', encoding='utf-8'
        )
        analysis = GOAnalysis(self.go_config())

        with mock.patch(
            'src.modules.go_analysis.subprocess.run',
            return_value=self.completed,
        ) as run_mock:
            analysis.add_go_annotations(str(input_path))
            analysis.compute_go_kfi(str(
                self.intermediate_dir / 'annotation_with_guides_go.tsv'
            ))
            analysis.count_go_parent_terms(str(fully_path))

        add_command, compute_command, parent_command = [
            call.args[0] for call in run_mock.call_args_list
        ]
        for command in (add_command, compute_command, parent_command):
            self.assert_workspace_command(command, '/DATA')

        self.assertEqual(add_command[-2:], [
            '/DATA/annotation.tsv',
            '/DATA/annotation_with_guides_go.tsv',
        ])
        self.assertEqual(
            compute_command[compute_command.index('-i') + 1],
            '/DATA/annotation_with_guides_go.tsv',
        )
        self.assertEqual(
            compute_command[compute_command.index('-o') + 1],
            '/DATA/go_kfi_table.csv',
        )
        self.assertEqual(
            compute_command[compute_command.index('-s') + 1],
            '/DATA/go_kfi_stat.csv',
        )
        self.assertEqual(
            parent_command[parent_command.index('--go_table') + 1],
            '/DATA/fully_targetable.csv',
        )
        self.assertEqual(
            parent_command[parent_command.index('--out_csv') + 1],
            '/DATA/go_parent_counts.csv',
        )
        self.assertEqual(
            parent_command[parent_command.index('--out_png') + 1],
            '/DATA/go_parent_counts.png',
        )

    def test_crisflash_uses_intermediate_workspace_for_all_io(self):
        config = {
            'output': {'intermediate': str(self.intermediate_dir)},
            'crispr': {
                'memory': '1g',
                'threads': 2,
                'docker_image': 'crisflash:test',
                'pam_pattern': 'NGG',
            },
        }
        genome_path = self.intermediate_dir / 'reference_for_offtarget.fa'
        gene_path = self.intermediate_dir / 'single_genes.fa'
        gene_path.write_text('>Single_1\nACGT\n', encoding='utf-8')
        (self.intermediate_dir / 'single_genes_targets.bed').touch()
        completed = subprocess.CompletedProcess(
            args=[],
            returncode=0,
            stdout=(
                "[crisflash] 0 gRNA sequences discovered in 'Single_1'."
            ),
            stderr='',
        )

        with mock.patch(
            'src.modules.guide_design.subprocess.run',
            return_value=completed,
        ) as run_mock:
            GuideDesign(config).run_crisflash(
                str(genome_path), str(gene_path), 'single_genes'
            )

        command = run_mock.call_args.args[0]
        self.assert_workspace_command(command, '/data')
        self.assertEqual(
            command[command.index('-g') + 1],
            '/data/reference_for_offtarget.fa',
        )
        self.assertEqual(
            command[command.index('-s') + 1],
            '/data/single_genes.fa',
        )
        self.assertEqual(
            command[command.index('-o') + 1],
            '/data/single_genes_targets.bed',
        )
        crisflash_arguments = command[command.index('crisflash') + 1:]
        self.assertEqual(
            crisflash_arguments[crisflash_arguments.index('-m') + 1],
            '2',
        )

    def test_seqkit_uses_intermediate_workspace_for_absolute_protein_path(self):
        config = {
            'output': {'intermediate': str(self.intermediate_dir)},
        }
        protein_path = self.intermediate_dir / 'query.pep'

        with mock.patch(
            'src.modules.gene_prediction.subprocess.run',
            return_value=self.completed,
        ) as run_mock:
            result = GenePrediction(config).get_sequence_stats(
                str(protein_path)
            )

        command = run_mock.call_args.args[0]
        self.assert_workspace_command(command, '/seqkit_work')
        self.assertEqual(
            command[command.index('-w') + 1], '/seqkit_work'
        )
        self.assertEqual(command[-1], '/seqkit_work/query.pep')
        self.assertEqual(result, 'completed')


if __name__ == '__main__':
    unittest.main()
