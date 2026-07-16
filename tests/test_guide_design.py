import csv
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from src.modules.guide_design import GuideDesign


class SiteSpecificTargetTests(unittest.TestCase):
    def test_run_publishes_normalized_generation_without_per_type_join(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            gene_path = intermediate / 'single_genes.fa'
            target_path = intermediate / 'single_targets.bed'
            specific_path = intermediate / 'single_specific.bed'
            annotation_path = intermediate / 'annotation_classified.txt'
            annotation_path.write_text('header\n', encoding='utf-8')
            published = {
                'normalized_genes': str(intermediate / 'generation/genes.tsv'),
                'normalized_guides': str(intermediate / 'generation/guides.tsv'),
                'normalized_gene_guide_manifest': str(
                    intermediate / 'generation/manifest.json'
                ),
                'annotation_with_guides_final': str(
                    intermediate / 'generation/legacy.tsv'
                ),
                'normalized_gene_guide_current': str(
                    intermediate / 'current.json'
                ),
                'generation_dir': str(intermediate / 'generation'),
            }
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })
            guide_design.prepare_genome_assembly = mock.Mock(
                return_value=str(intermediate / 'reference.fa')
            )
            guide_design.prepare_gene_sequences = mock.Mock(
                return_value={'single_genes': str(gene_path)}
            )
            guide_design.run_crisflash = mock.Mock(
                return_value=str(target_path)
            )
            guide_design.filter_site_specific_targets = mock.Mock(
                return_value=str(specific_path)
            )
            guide_design.merge_with_annotation = mock.Mock(
                side_effect=AssertionError('legacy per-type join was called')
            )
            guide_design.create_final_combined_table = mock.Mock(
                side_effect=AssertionError('legacy combined copy was called')
            )

            with mock.patch(
                'src.modules.guide_design.'
                'publish_normalized_gene_guide_generation',
                return_value=published,
            ) as publish_mock:
                results = guide_design.run()

            publish_mock.assert_called_once_with(
                annotation_path=str(annotation_path),
                output_dir=temp_dir,
                multi_bed_path=None,
                single_bed_path=str(specific_path),
            )
            guide_design.merge_with_annotation.assert_not_called()
            guide_design.create_final_combined_table.assert_not_called()
            self.assertEqual(
                results['normalized_genes'], published['normalized_genes']
            )
            self.assertEqual(
                results['normalized_guides'], published['normalized_guides']
            )
            self.assertEqual(
                results['final_table'],
                published['annotation_with_guides_final'],
            )
            self.assertEqual(
                results['target_results']['single_genes'][
                    'annotation_with_guides'
                ],
                published['annotation_with_guides_final'],
            )

    def test_empty_alignment_directory_rejects_expected_multi_genes(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'multi_genes.fa').write_text(
                '>Multi_1_tx1\nACGT\n>Multi_1_tx2\nACGA\n',
                encoding='utf-8',
            )
            (intermediate / 'multi_genes_aligned').mkdir()
            (intermediate / 'single_genes.fa').write_text(
                '>Single_1\nACGT\n', encoding='utf-8'
            )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            with self.assertRaisesRegex(
                RuntimeError, 'missing=Multi_1'
            ):
                guide_design.prepare_gene_sequences()

    def test_partial_alignment_output_rejects_missing_multi_group(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'multi_genes.fa').write_text(
                '>Multi_1_tx1\nACGT\n>Multi_2_tx1\nTGCA\n',
                encoding='utf-8',
            )
            (intermediate / 'annotation_classified.txt').write_text(
                'classification\tunique_id\n'
                'Multi\tMulti_1\nMulti\tMulti_2\n',
                encoding='utf-8',
            )
            aligned_dir = intermediate / 'multi_genes_aligned'
            aligned_dir.mkdir()
            (aligned_dir / 'Multi_1_aligned_consensus.fa').write_text(
                '>Multi_1_aligned_consensus\nACGT\n', encoding='utf-8'
            )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            with self.assertRaisesRegex(
                RuntimeError, 'missing=Multi_2'
            ):
                guide_design.prepare_gene_sequences()

    def test_complete_alignment_output_prepares_all_multi_groups(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'multi_genes.fa').write_text(
                '>Multi_1_tx1\nACGT\n>Multi_2_tx1\nTGCA\n',
                encoding='utf-8',
            )
            aligned_dir = intermediate / 'multi_genes_aligned'
            aligned_dir.mkdir()
            for multi_id, sequence in (
                ('Multi_1', 'ACGT'), ('Multi_2', 'TGCA')
            ):
                (aligned_dir / f'{multi_id}_aligned_consensus.fa').write_text(
                    f'>{multi_id} consensus sequence\n{sequence}\n',
                    encoding='utf-8',
                )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            prepared_files = guide_design.prepare_gene_sequences()

            combined_path = Path(prepared_files['multi_genes'])
            self.assertEqual(
                combined_path.read_text(encoding='utf-8'),
                '>Multi_1_aligned_consensus\nACGT\n'
                '>Multi_2_aligned_consensus\nTGCA\n',
            )

    def test_empty_consensus_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'multi_genes.fa').write_text(
                '>Multi_1_tx1\nACGT\n', encoding='utf-8'
            )
            aligned_dir = intermediate / 'multi_genes_aligned'
            aligned_dir.mkdir()
            (aligned_dir / 'Multi_1_aligned_consensus.fa').touch()
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            with self.assertRaisesRegex(
                RuntimeError, 'Invalid aligned consensus FASTA'
            ):
                guide_design.prepare_gene_sequences()

    def test_declared_multi_group_without_extracted_sequence_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'annotation_classified.txt').write_text(
                'classification\tunique_id\nMulti\tMulti_1\n',
                encoding='utf-8',
            )
            (intermediate / 'multi_genes.fa').touch()
            (intermediate / 'multi_genes_aligned').mkdir()
            (intermediate / 'single_genes.fa').write_text(
                '>Single_1\nACGT\n', encoding='utf-8'
            )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            with self.assertRaisesRegex(
                RuntimeError, 'missing=Multi_1'
            ):
                guide_design.prepare_gene_sequences()

    def test_empty_alignment_directory_is_allowed_when_multi_count_is_zero(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'multi_genes_aligned').mkdir()
            (intermediate / 'single_genes.fa').write_text(
                '>Single_1\nACGT\n', encoding='utf-8'
            )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            prepared_files = guide_design.prepare_gene_sequences()

            self.assertEqual(set(prepared_files), {'single_genes'})

    def test_declared_single_gene_without_extracted_sequence_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'annotation_classified.txt').write_text(
                'classification\tunique_id\nSingle\tSingle_1\n',
                encoding='utf-8',
            )
            (intermediate / 'single_genes.fa').touch()
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            with self.assertRaisesRegex(RuntimeError, 'missing=Single_1'):
                guide_design.prepare_gene_sequences()

    def test_partial_single_gene_fasta_rejects_missing_declaration(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'annotation_classified.txt').write_text(
                'classification\tunique_id\n'
                'Single\tSingle_1\nSingle\tSingle_2\n',
                encoding='utf-8',
            )
            (intermediate / 'single_genes.fa').write_text(
                '>Single_1\nACGT\n', encoding='utf-8'
            )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            with self.assertRaisesRegex(RuntimeError, 'missing=Single_2'):
                guide_design.prepare_gene_sequences()

    def test_stale_single_gene_fasta_record_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'annotation_classified.txt').write_text(
                'classification\tunique_id\nSingle\tSingle_1\n',
                encoding='utf-8',
            )
            (intermediate / 'single_genes.fa').write_text(
                '>Single_1\nACGT\n>Single_stale\nTGCA\n',
                encoding='utf-8',
            )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            with self.assertRaisesRegex(
                RuntimeError, 'unexpected=Single_stale'
            ):
                guide_design.prepare_gene_sequences()

    def test_empty_single_gene_sequence_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'annotation_classified.txt').write_text(
                'classification\tunique_id\nSingle\tSingle_1\n',
                encoding='utf-8',
            )
            (intermediate / 'single_genes.fa').write_text(
                '>Single_1\n', encoding='utf-8'
            )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            with self.assertRaisesRegex(RuntimeError, 'empty sequence'):
                guide_design.prepare_gene_sequences()

    def test_matching_single_gene_declarations_are_prepared(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'annotation_classified.txt').write_text(
                'classification\tunique_id\n'
                'Single\tSingle_1\nSingle\tSingle_2\n',
                encoding='utf-8',
            )
            (intermediate / 'single_genes.fa').write_text(
                '>Single_1\nACGT\n>Single_2\nTGCA\n',
                encoding='utf-8',
            )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            prepared_files = guide_design.prepare_gene_sequences()

            self.assertEqual(set(prepared_files), {'single_genes'})

    def test_missing_optional_genome_assembly_uses_transcriptome(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            transcriptome_path = temp_path / 'transcriptome.fa'
            transcriptome_path.write_text(
                '>tx1\nACGTRY\n', encoding='utf-8'
            )
            guide_design = GuideDesign({
                'input': {'transcriptome': str(transcriptome_path)},
                'output': {'intermediate': temp_dir},
            })

            output_path = guide_design.prepare_genome_assembly()

            self.assertEqual(
                Path(output_path).read_text(encoding='utf-8'),
                '>tx1\nACGTNN\n',
            )

    def test_site_specific_target_requires_exact_mismatch_counts(self):
        cases = {
            'ACGTNGG/1/0': True,
            'ACGTNGG/1/45': False,
            'ACGTNGG/0/0': False,
            'ACGTNGG/2/0': False,
            'ACGTNGG/11/0': False,
            'ACGTNGG/1/not-a-number': False,
            'ACGTNGG/1/0/extra': False,
            '/1/0': False,
            'ACGTNGG': False,
        }

        for target_field, expected in cases.items():
            with self.subTest(target_field=target_field):
                self.assertEqual(
                    GuideDesign._is_site_specific_target(target_field),
                    expected,
                )

    def test_filter_keeps_only_mm0_one_and_mm1_2_zero(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / 'multi_genes_targets.bed'
            input_path.write_text(
                ''.join([
                    'Multi_1\t10\t33\tACGTNGG/1/0\t1.0\t+\n',
                    'Multi_1\t20\t43\tACGTNGG/1/45\t0.9\t+\n',
                    'Multi_1\t30\t53\tACGTNGG/0/0\t0.8\t-\n',
                    'Multi_1\t40\t63\tACGTNGG/2/0\t0.7\t-\n',
                    'Multi_1\t50\t73\tACGTNGG/1/bad\t0.6\t+\n',
                ]),
                encoding='utf-8',
            )

            guide_design = GuideDesign({
                'output': {
                    'intermediate': temp_dir,
                },
            })
            output_path = guide_design.filter_site_specific_targets(
                str(input_path),
                'multi_genes',
            )

            self.assertEqual(
                Path(output_path).read_text(encoding='utf-8'),
                'Multi_1\t10\t33\tACGTNGG/1/0\t1.0\t+\n',
            )

    def test_merge_preserves_empty_identifier_columns_before_guides(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            intermediate = Path(temp_dir)
            (intermediate / 'annotation_classified.txt').write_text(
                'query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\t'
                'description\tgene_name\tclassification\tunique_id\t'
                'annotation_source_database\tsgd_id\tensembl_id\n'
                'q1\tSGDID:S000000056\t\t383\t0\t154.8\tdesc\tBDH1\t'
                'Single\tSingle_0000000001\tSGD\tS000000056\t\n',
                encoding='utf-8',
            )
            targets_path = intermediate / 'single_genes_targets.bed'
            targets_path.write_text(
                'Single_0000000001\t1\t23\tACGT\t0\t0\n',
                encoding='utf-8',
            )
            guide_design = GuideDesign({
                'output': {'intermediate': temp_dir},
            })

            output_path = guide_design.merge_with_annotation(
                str(targets_path), 'single_genes'
            )
            with open(output_path, encoding='utf-8', newline='') as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['ensembl_id'], '')
            self.assertEqual(
                rows[0]['guide_targets'],
                'Single_0000000001/1/23/ACGT/0/0',
            )


if __name__ == '__main__':
    unittest.main()
