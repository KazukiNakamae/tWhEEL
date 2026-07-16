import tempfile
import unittest
from pathlib import Path
from unittest import mock

from src.modules.annotation import Annotation
from src.modules.sgd_annotation import Sgd_annotation
from src.modules.transcript_mapping import TranscriptMappingError


ANNOTATION_CLASSES = (Annotation, Sgd_annotation)


class TranscriptMappingContractTests(unittest.TestCase):
    @staticmethod
    def _config(temp_dir, mapping_path):
        return {
            'input': {
                'protein_transcript_mapping': str(mapping_path),
            },
            'output': {
                'intermediate': str(temp_dir),
            },
        }

    @staticmethod
    def _write_inputs(
        temp_dir,
        mapping_content,
        protein_content='>protein1\nMPEPTIDE\n',
        transcript_content='>tx1\nATGC\n',
    ):
        mapping_path = temp_dir / 'mapping.tsv'
        protein_path = temp_dir / 'query.pep'
        transcript_path = temp_dir / 'transcriptome.fa'
        mapping_path.write_text(mapping_content, encoding='utf-8')
        protein_path.write_text(protein_content, encoding='utf-8')
        transcript_path.write_text(transcript_content, encoding='utf-8')
        return mapping_path, protein_path, transcript_path

    def test_valid_mapping_is_loaded_and_validated_by_both_classes(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            mapping_path, protein_path, transcript_path = self._write_inputs(
                temp_dir,
                '# transcript-to-protein mapping\n'
                'tx1\tprotein1\n'
                'tx2\tprotein2\n',
                protein_content=(
                    '>protein1 optional description\nMPEPTIDE\n'
                    '>protein2\nMSECOND\n'
                ),
                transcript_content=(
                    '>tx1 description\nATGC\n>tx2\nGCTA\n'
                    '>unmapped_noncoding_transcript\nAAAA\n'
                ),
            )

            for annotation_class in ANNOTATION_CLASSES:
                with self.subTest(annotation_class=annotation_class.__name__):
                    annotation = annotation_class(
                        self._config(temp_dir, mapping_path)
                    )

                    mapping = annotation.load_protein_transcript_mapping()
                    validation = annotation.validate_protein_transcript_mapping(
                        protein_path,
                        transcript_path,
                    )

                    self.assertEqual(
                        mapping,
                        {'protein1': 'tx1', 'protein2': 'tx2'},
                    )
                    self.assertEqual(validation['protein_count'], 2)
                    self.assertEqual(validation['mapped_transcript_count'], 2)
                    self.assertEqual(validation['transcript_count'], 3)

    def test_run_rejects_missing_transcript_before_fasta_search(self):
        for annotation_class in ANNOTATION_CLASSES:
            with self.subTest(annotation_class=annotation_class.__name__):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    temp_dir = Path(temp_dir_name)
                    mapping_path, protein_path, _ = self._write_inputs(
                        temp_dir,
                        'tx_missing\tprotein1\n',
                    )
                    config = self._config(temp_dir, mapping_path)
                    config['input']['protein_sequences'] = str(protein_path)
                    config['transdecoder'] = {
                        'skip_if_protein_provided': True,
                    }
                    annotation = annotation_class(config)

                    with mock.patch.object(
                        annotation,
                        'run_fasta_search',
                    ) as search_mock:
                        with self.assertRaisesRegex(
                            TranscriptMappingError,
                            'protein1->tx_missing',
                        ):
                            annotation.run()

                    search_mock.assert_not_called()

    def test_missing_query_protein_mapping_is_rejected_by_both_classes(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            mapping_path, protein_path, transcript_path = self._write_inputs(
                temp_dir,
                'tx1\tprotein1\n',
                protein_content='>protein1\nMPEPTIDE\n>protein2\nMSECOND\n',
                transcript_content='>tx1\nATGC\n>tx2\nGCTA\n',
            )

            for annotation_class in ANNOTATION_CLASSES:
                with self.subTest(annotation_class=annotation_class.__name__):
                    annotation = annotation_class(
                        self._config(temp_dir, mapping_path)
                    )
                    annotation.load_protein_transcript_mapping()

                    with self.assertRaises(TranscriptMappingError) as context:
                        annotation.validate_protein_transcript_mapping(
                            protein_path,
                            transcript_path,
                        )

                    message = str(context.exception)
                    self.assertIn('protein ID(s) without a mapping', message)
                    self.assertIn('protein2', message)
                    self.assertIn('rerun the annotation step', message)

    def test_mapping_protein_absent_from_query_is_rejected_by_both_classes(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            mapping_path, protein_path, transcript_path = self._write_inputs(
                temp_dir,
                'tx1\tprotein1\ntx2\tprotein2\n',
                transcript_content='>tx1\nATGC\n>tx2\nGCTA\n',
            )

            for annotation_class in ANNOTATION_CLASSES:
                with self.subTest(annotation_class=annotation_class.__name__):
                    annotation = annotation_class(
                        self._config(temp_dir, mapping_path)
                    )
                    annotation.load_protein_transcript_mapping()

                    with self.assertRaises(TranscriptMappingError) as context:
                        annotation.validate_protein_transcript_mapping(
                            protein_path,
                            transcript_path,
                        )

                    message = str(context.exception)
                    self.assertIn(
                        'mapping protein ID(s) absent from query.pep',
                        message,
                    )
                    self.assertIn('protein2', message)

    def test_mapped_transcript_absent_from_fasta_is_rejected_by_both_classes(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            mapping_path, protein_path, transcript_path = self._write_inputs(
                temp_dir,
                'tx_missing\tprotein1\n',
            )

            for annotation_class in ANNOTATION_CLASSES:
                with self.subTest(annotation_class=annotation_class.__name__):
                    annotation = annotation_class(
                        self._config(temp_dir, mapping_path)
                    )
                    annotation.load_protein_transcript_mapping()

                    with self.assertRaises(TranscriptMappingError) as context:
                        annotation.validate_protein_transcript_mapping(
                            protein_path,
                            transcript_path,
                        )

                    message = str(context.exception)
                    self.assertIn(
                        'mapped transcript ID(s) absent from transcriptome.fa',
                        message,
                    )
                    self.assertIn('protein1->tx_missing', message)

    def test_malformed_empty_and_headered_mappings_are_rejected(self):
        cases = (
            (
                'malformed',
                'tx1 protein1\n',
                'expected a headerless transcript_id<TAB>protein_id row',
            ),
            (
                'empty',
                '\n# comments are not data rows\n',
                'contains no data rows',
            ),
            (
                'header',
                'transcript_id\tprotein_id\ntx1\tprotein1\n',
                'header rows are not allowed',
            ),
        )

        for name, mapping_content, expected_message in cases:
            for annotation_class in ANNOTATION_CLASSES:
                with self.subTest(
                    case=name,
                    annotation_class=annotation_class.__name__,
                ):
                    with tempfile.TemporaryDirectory() as temp_dir_name:
                        temp_dir = Path(temp_dir_name)
                        mapping_path, _, _ = self._write_inputs(
                            temp_dir,
                            mapping_content,
                        )
                        annotation = annotation_class(
                            self._config(temp_dir, mapping_path)
                        )

                        with self.assertRaisesRegex(
                            TranscriptMappingError,
                            expected_message,
                        ):
                            annotation.load_protein_transcript_mapping()

    def test_malformed_row_after_ten_valid_rows_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            mapping_path, _, _ = self._write_inputs(
                temp_dir,
                ('tx1\tprotein1\n' * 10) + 'malformed_row\n',
            )

            with self.assertRaisesRegex(
                TranscriptMappingError,
                r':11: expected a headerless',
            ):
                Annotation(
                    self._config(temp_dir, mapping_path)
                ).load_protein_transcript_mapping()

    def test_conflicting_duplicate_mapping_is_rejected_by_both_classes(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            mapping_path, _, _ = self._write_inputs(
                temp_dir,
                'tx1\tprotein1\ntx2\tprotein1\n',
            )

            for annotation_class in ANNOTATION_CLASSES:
                with self.subTest(annotation_class=annotation_class.__name__):
                    annotation = annotation_class(
                        self._config(temp_dir, mapping_path)
                    )

                    with self.assertRaises(TranscriptMappingError) as context:
                        annotation.load_protein_transcript_mapping()

                    message = str(context.exception)
                    self.assertIn(
                        'Conflicting protein-transcript mappings',
                        message,
                    )
                    self.assertIn("protein_id 'protein1'", message)
                    self.assertIn("'tx1' at line 1", message)
                    self.assertIn("'tx2' at line 2", message)

    def test_identical_duplicate_mapping_is_allowed_by_both_classes(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            mapping_path, protein_path, transcript_path = self._write_inputs(
                temp_dir,
                'tx1\tprotein1\ntx1\tprotein1\n',
            )

            for annotation_class in ANNOTATION_CLASSES:
                with self.subTest(annotation_class=annotation_class.__name__):
                    annotation = annotation_class(
                        self._config(temp_dir, mapping_path)
                    )

                    self.assertEqual(
                        annotation.load_protein_transcript_mapping(),
                        {'protein1': 'tx1'},
                    )
                    annotation.validate_protein_transcript_mapping(
                        protein_path,
                        transcript_path,
                    )

    def test_annotation_query_absent_from_mapping_is_not_treated_as_untargetable(self):
        annotation_table = (
            'query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\t'
            'description\n'
            'protein_missing\tsubject1\t\t8\t0\t10\tExample protein\n'
        )

        for annotation_class in ANNOTATION_CLASSES:
            with self.subTest(annotation_class=annotation_class.__name__):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    temp_dir = Path(temp_dir_name)
                    annotation_path = temp_dir / 'annotation_table.txt'
                    classified_path = temp_dir / 'annotation_classified.txt'
                    annotation_path.write_text(
                        annotation_table,
                        encoding='utf-8',
                    )
                    classified_path.write_text(
                        'existing classified output\n',
                        encoding='utf-8',
                    )
                    annotation = annotation_class({
                        'output': {'intermediate': str(temp_dir)},
                    })
                    annotation.protein_transcript_mapping = {
                        'protein1': 'tx1',
                    }

                    with self.assertRaisesRegex(
                        TranscriptMappingError,
                        'Refusing to classify this entry as untargetable',
                    ):
                        annotation.classify_genes(annotation_path)

                    self.assertEqual(
                        classified_path.read_text(encoding='utf-8'),
                        'existing classified output\n',
                    )

    def test_malformed_annotation_row_is_rejected_in_mapping_mode(self):
        annotation_table = (
            'query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\t'
            'description\nprotein1\n'
        )

        for annotation_class in ANNOTATION_CLASSES:
            with self.subTest(annotation_class=annotation_class.__name__):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    temp_dir = Path(temp_dir_name)
                    annotation_path = temp_dir / 'annotation_table.txt'
                    classified_path = temp_dir / 'annotation_classified.txt'
                    annotation_path.write_text(
                        annotation_table,
                        encoding='utf-8',
                    )
                    classified_path.write_text(
                        'existing classified output\n',
                        encoding='utf-8',
                    )
                    annotation = annotation_class({
                        'output': {'intermediate': str(temp_dir)},
                    })
                    annotation.protein_transcript_mapping = {
                        'protein1': 'tx1',
                    }

                    with self.assertRaisesRegex(
                        TranscriptMappingError,
                        'Malformed annotation row',
                    ):
                        annotation.classify_genes(annotation_path)

                    self.assertEqual(
                        classified_path.read_text(encoding='utf-8'),
                        'existing classified output\n',
                    )

    def test_missing_transcript_during_extraction_preserves_existing_output(self):
        for annotation_class in ANNOTATION_CLASSES:
            with self.subTest(annotation_class=annotation_class.__name__):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    temp_dir = Path(temp_dir_name)
                    transcript_path = temp_dir / 'transcriptome.fa'
                    output_path = temp_dir / 'multi_genes.fa'
                    transcript_path.write_text(
                        '>tx1\nATGC\n',
                        encoding='utf-8',
                    )
                    output_path.write_text(
                        'existing complete output\n',
                        encoding='utf-8',
                    )
                    annotation = annotation_class({
                        'output': {'intermediate': str(temp_dir)},
                    })
                    annotation.protein_transcript_mapping = {
                        'protein1': 'tx1',
                        'protein2': 'tx_missing',
                    }

                    with self.assertRaises(TranscriptMappingError) as context:
                        annotation._extract_sequences_by_ids(
                            transcript_path,
                            output_path,
                            {
                                'Multi_0000000001': [
                                    'protein1',
                                    'protein2',
                                ],
                            },
                        )

                    self.assertIn(
                        'protein2->tx_missing',
                        str(context.exception),
                    )
                    self.assertEqual(
                        output_path.read_text(encoding='utf-8'),
                        'existing complete output\n',
                    )

    def test_all_gene_sequences_are_resolved_before_any_output_is_created(self):
        classified_table = (
            'query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\t'
            'description\tgene_name\tclassification\tunique_id\n'
            'protein1\tsubject1\t\t8\t0\t10\tExample\tGENE1\tMulti\t'
            'Multi_0000000001\n'
            'protein2\tsubject2\t\t8\t0\t10\tExample\tGENE2\tSingle\t'
            'Single_0000000001\n'
        )

        for annotation_class in ANNOTATION_CLASSES:
            with self.subTest(annotation_class=annotation_class.__name__):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    temp_dir = Path(temp_dir_name)
                    classified_path = temp_dir / 'annotation_classified.txt'
                    classified_path.write_text(
                        classified_table,
                        encoding='utf-8',
                    )
                    (temp_dir / 'transcriptome.fa').write_text(
                        '>tx1\nATGC\n',
                        encoding='utf-8',
                    )
                    annotation = annotation_class({
                        'output': {'intermediate': str(temp_dir)},
                    })
                    annotation.protein_transcript_mapping = {
                        'protein1': 'tx1',
                        'protein2': 'tx_missing',
                    }

                    with self.assertRaisesRegex(
                        TranscriptMappingError,
                        'protein2->tx_missing',
                    ):
                        annotation.extract_gene_sequences(classified_path)

                    self.assertFalse((temp_dir / 'multi_genes.fa').exists())
                    self.assertFalse((temp_dir / 'single_genes.fa').exists())
                    self.assertFalse((temp_dir / 'multi_genes').exists())

    def test_transdecoder_resolution_removes_only_a_terminal_orf_suffix(self):
        for annotation_class in ANNOTATION_CLASSES:
            with self.subTest(annotation_class=annotation_class.__name__):
                annotation = annotation_class({
                    'output': {'intermediate': '/unused'},
                })

                self.assertEqual(
                    annotation._resolve_sequence_id('sample.part.p2'),
                    'sample.part',
                )
                self.assertEqual(
                    annotation._resolve_sequence_id('sample.part'),
                    'sample.part',
                )

                with self.assertRaisesRegex(
                    TranscriptMappingError,
                    'TransDecoder-style',
                ):
                    annotation._resolve_required_sequences(
                        {'Single_0000000001': 'sample.part.p2'},
                        {},
                    )

    def test_mapping_mode_never_falls_back_to_transdecoder_id_rules(self):
        mapping_states = (
            ('empty mapping state', {}),
            ('unrelated mapping entry', {'other_protein': 'other_tx'}),
        )

        for annotation_class in ANNOTATION_CLASSES:
            for state_name, mapping in mapping_states:
                with self.subTest(
                    annotation_class=annotation_class.__name__,
                    state=state_name,
                ):
                    annotation = annotation_class({
                        'output': {'intermediate': '/unused'},
                    })
                    annotation.protein_transcript_mapping = mapping

                    with self.assertRaisesRegex(
                        TranscriptMappingError,
                        'Protein query ID is absent from the loaded',
                    ):
                        annotation._resolve_sequence_id('tx1.p1')


if __name__ == '__main__':
    unittest.main()
