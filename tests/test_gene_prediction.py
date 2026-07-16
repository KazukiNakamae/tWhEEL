import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src.modules.gene_prediction import GenePrediction


class TransDecoderMinimumProteinLengthTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir_context = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temp_dir_context.name)
        self.intermediate_dir = self.temp_dir / 'intermediate'
        self.intermediate_dir.mkdir()
        self.input_path = self.temp_dir / 'transcriptome.fa'
        self.input_path.write_text('>tx1\nATGGCC\n', encoding='utf-8')
        self.config = {
            'output': {'intermediate': str(self.intermediate_dir)},
            'transdecoder': {
                'docker_image': 'transdecoder:test',
                'min_protein_len': 120,
            },
        }

    def tearDown(self):
        self.temp_dir_context.cleanup()

    def test_longorfs_receives_configured_minimum_protein_length(self):
        commands = []

        def run_command(command, **kwargs):
            commands.append(command)
            if 'TransDecoder.Predict' in command:
                output_path = (
                    self.intermediate_dir
                    / 'transdecoder_output'
                    / 'transcriptome.fa.transdecoder.pep'
                )
                output_path.write_text('>tx1.p1\nMA\n', encoding='utf-8')
            return subprocess.CompletedProcess(command, 0, '', '')

        with patch(
            'src.modules.gene_prediction.subprocess.run',
            side_effect=run_command,
        ):
            result_path = GenePrediction(self.config).run_transdecoder(
                str(self.input_path)
            )

        longorfs_command = next(
            command for command in commands
            if 'TransDecoder.LongOrfs' in command
        )
        predict_command = next(
            command for command in commands
            if 'TransDecoder.Predict' in command
        )
        minimum_index = longorfs_command.index('-m')
        self.assertEqual(longorfs_command[minimum_index + 1], '120')
        self.assertNotIn('-m', predict_command)
        self.assertTrue(Path(result_path).is_file())

    def test_missing_setting_uses_transdecoder_default(self):
        del self.config['transdecoder']['min_protein_len']

        minimum = GenePrediction(
            self.config
        )._get_transdecoder_min_protein_len()

        self.assertEqual(minimum, 100)

    def test_invalid_minimum_protein_length_is_rejected(self):
        for value in (0, -1, True, 85.5, '', 'not-an-integer'):
            with self.subTest(value=value):
                self.config['transdecoder']['min_protein_len'] = value

                with self.assertRaisesRegex(
                    ValueError, 'must be a positive integer'
                ):
                    GenePrediction(
                        self.config
                    )._get_transdecoder_min_protein_len()


if __name__ == '__main__':
    unittest.main()
