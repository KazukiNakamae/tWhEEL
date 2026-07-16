import tempfile
import unittest
from pathlib import Path

from src.modules.annotation import Annotation
from src.modules.gene_prediction import GenePrediction
from src.modules.input_mode import (
    EXISTING_PROTEIN_MODE,
    TRANSDECODER_MODE,
    determine_protein_input_mode,
)
from src.modules.sgd_annotation import Sgd_annotation


class ProteinInputModeTests(unittest.TestCase):
    @staticmethod
    def _config(temp_dir, transcriptome_path, protein_path, mapping_path, skip):
        return {
            'input': {
                'transcriptome': str(transcriptome_path),
                'protein_sequences': protein_path,
                'protein_transcript_mapping': mapping_path,
            },
            'output': {
                'intermediate': str(temp_dir),
            },
            'transdecoder': {
                'skip_if_protein_provided': skip,
            },
        }

    def test_all_modules_resolve_the_same_effective_mode(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            transcriptome_path = temp_dir / 'transcriptome.fasta'
            protein_path = temp_dir / 'proteins.fasta'
            mapping_path = temp_dir / 'mapping.tsv'
            missing_protein_path = temp_dir / 'missing.fasta'
            transcriptome_path.write_text('>tx1\nATGC\n', encoding='utf-8')
            protein_path.write_text('>protein1\nMPEPTIDE\n', encoding='utf-8')
            mapping_path.write_text('protein1\ttx1\n', encoding='utf-8')

            cases = [
                ('', True, '', TRANSDECODER_MODE),
                (str(missing_protein_path), True, '', TRANSDECODER_MODE),
                (str(protein_path), False, '', TRANSDECODER_MODE),
                (
                    str(protein_path),
                    True,
                    str(mapping_path),
                    EXISTING_PROTEIN_MODE,
                ),
            ]

            for configured_protein, skip, configured_mapping, expected in cases:
                with self.subTest(
                    protein=configured_protein,
                    skip=skip,
                    expected=expected,
                ):
                    config = self._config(
                        temp_dir,
                        transcriptome_path,
                        configured_protein,
                        configured_mapping,
                        skip,
                    )
                    annotation = Annotation(config)
                    sgd_annotation = Sgd_annotation(config)

                    self.assertEqual(
                        determine_protein_input_mode(config),
                        expected,
                    )
                    self.assertEqual(
                        GenePrediction(config).determine_input_mode(),
                        expected,
                    )
                    self.assertEqual(
                        annotation._validate_existing_protein_mode_requirements(),
                        expected,
                    )
                    self.assertEqual(
                        sgd_annotation._validate_existing_protein_mode_requirements(),
                        expected,
                    )

    def test_mapping_is_required_only_for_effective_existing_protein_mode(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            transcriptome_path = temp_dir / 'transcriptome.fasta'
            protein_path = temp_dir / 'proteins.fasta'
            transcriptome_path.write_text('>tx1\nATGC\n', encoding='utf-8')
            protein_path.write_text('>protein1\nMPEPTIDE\n', encoding='utf-8')
            config = self._config(
                temp_dir,
                transcriptome_path,
                str(protein_path),
                '',
                True,
            )

            for annotation_class in (Annotation, Sgd_annotation):
                with self.subTest(annotation_class=annotation_class.__name__):
                    with self.assertRaisesRegex(
                        ValueError,
                        'protein_transcript_mapping',
                    ):
                        validator = annotation_class(
                            config
                        )._validate_existing_protein_mode_requirements
                        validator()

    def test_skip_flag_defaults_to_enabled(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            protein_path = temp_dir / 'proteins.fasta'
            protein_path.write_text('>protein1\nMPEPTIDE\n', encoding='utf-8')
            config = {
                'input': {
                    'protein_sequences': str(protein_path),
                },
            }

            self.assertEqual(
                determine_protein_input_mode(config),
                EXISTING_PROTEIN_MODE,
            )


if __name__ == '__main__':
    unittest.main()
