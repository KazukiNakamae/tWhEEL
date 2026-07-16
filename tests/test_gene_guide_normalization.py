import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from src.modules import gene_guide_normalization as normalization


ANNOTATION_COLUMNS = [
    'query_id', 'subject_id', 'identity', 'length', 'evalue', 'bit_score',
    'description', 'gene_name', 'classification', 'unique_id',
    'annotation_source_database', 'sgd_id', 'ensembl_id',
]


class GeneGuideNormalizationTests(unittest.TestCase):
    def setUp(self):
        self.temp_context = tempfile.TemporaryDirectory()
        self.temp_dir = Path(self.temp_context.name)
        self.annotation_path = self.temp_dir / 'annotation_classified.txt'

    def tearDown(self):
        self.temp_context.cleanup()

    @staticmethod
    def annotation_row(**overrides):
        row = {
            'query_id': 'query-single',
            'subject_id': 'sp|P00001|GENE_HUMAN',
            'identity': '90',
            'length': '100',
            'evalue': '1e-20',
            'bit_score': '120',
            'description': 'single description',
            'gene_name': 'GENE1',
            'classification': 'Single',
            'unique_id': 'Single_0000000001',
            'annotation_source_database': '',
            'sgd_id': '',
            'ensembl_id': '',
        }
        row.update(overrides)
        return row

    def write_annotation(self, rows):
        with self.annotation_path.open(
            'w', encoding='utf-8', newline=''
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=ANNOTATION_COLUMNS,
                delimiter='\t',
                lineterminator='\n',
            )
            writer.writeheader()
            writer.writerows(rows)

    def publish_single(self, bed_text=None):
        self.write_annotation([self.annotation_row()])
        bed_path = None
        if bed_text is not None:
            bed_path = self.temp_dir / 'single.bed'
            bed_path.write_text(bed_text, encoding='utf-8')
        return normalization.publish_normalized_gene_guide_generation(
            self.annotation_path,
            self.temp_dir,
            single_bed_path=bed_path,
        )

    def test_publishes_one_gene_row_and_one_typed_row_per_unique_guide(self):
        self.write_annotation([
            self.annotation_row(
                query_id='query-z',
                subject_id='sp|P00002|MULTI_HUMAN',
                classification='Multi',
                unique_id='Multi_0000000001',
                gene_name='MULTI',
                description='multi description',
                evalue='1e-10',
                bit_score='100',
            ),
            self.annotation_row(
                query_id='query-a',
                subject_id='sp|P00002|MULTI_HUMAN',
                classification='Multi',
                unique_id='Multi_0000000001',
                gene_name='MULTI',
                description='multi description',
                evalue='1e-50',
                bit_score='200',
            ),
            self.annotation_row(),
        ])
        multi_bed = self.temp_dir / 'multi.bed'
        multi_bed.write_text(
            'Multi_0000000001_aligned_consensus\t10\t33\tACGT/1/0\t1\t+\n'
            'Multi_0000000001_aligned_consensus\t10\t33\tACGT/1/0\t1.0\t+\n',
            encoding='utf-8',
        )
        single_bed = self.temp_dir / 'single.bed'
        single_bed.write_text(
            'Single_0000000001\t20\t43\tTGCA/1/0\t0.75\t-\n',
            encoding='utf-8',
        )

        published = normalization.publish_normalized_gene_guide_generation(
            self.annotation_path,
            self.temp_dir,
            multi_bed_path=multi_bed,
            single_bed_path=single_bed,
        )
        loaded = normalization.load_current_normalized_gene_guide_data(
            self.temp_dir
        )

        self.assertEqual(len(loaded['genes']), 2)
        self.assertEqual(len(loaded['guides']), 2)
        multi_gene = next(
            gene for gene in loaded['genes']
            if gene['classification'] == 'Multi'
        )
        self.assertEqual(multi_gene['query_ids'], ['query-a', 'query-z'])
        self.assertEqual(multi_gene['representative_query_id'], 'query-a')
        self.assertEqual(multi_gene['annotation_row_count'], 2)
        self.assertEqual(multi_gene['guide_count'], 1)
        multi_guide = next(
            guide for guide in loaded['guides']
            if guide['unique_id'] == multi_gene['unique_id']
        )
        self.assertEqual(
            multi_guide['reference'],
            'Multi_0000000001_aligned_consensus',
        )
        self.assertEqual(multi_guide['mm0'], 1)
        self.assertEqual(multi_guide['mm1_2'], 0)
        self.assertEqual(multi_guide['off_target_score'], 1.0)

        with Path(published['annotation_with_guides_final']).open(
            encoding='utf-8', newline=''
        ) as handle:
            legacy_rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(legacy_rows), 2)
        self.assertEqual(
            sum(
                row['unique_id'] == 'Multi_0000000001'
                for row in legacy_rows
            ),
            1,
        )
        self.assertFalse(
            (self.temp_dir / normalization.LEGACY_FINAL_FILENAME).exists()
        )
        self.assertTrue(
            (self.temp_dir / normalization.CURRENT_POINTER_FILENAME).is_file()
        )

    def test_parses_six_seven_and_eight_part_legacy_tokens(self):
        published = self.publish_single(
            'Single_0000000001/1/23/AAAA/1/0\n'
            'Single_0000000001/2/24/CCCC/1/0/0.5\n'
            'Single_0000000001/3/25/GGGG/1/0/0.25/-\n'
        )

        loaded = normalization.load_normalized_gene_guide_data(
            published['normalized_genes'],
            published['normalized_guides'],
            published['normalized_gene_guide_manifest'],
        )

        self.assertEqual(
            [guide['parse_status'] for guide in loaded['guides']],
            ['parsed_6_part', 'parsed_7_part', 'parsed_8_part'],
        )
        self.assertEqual(
            [guide['off_target_score'] for guide in loaded['guides']],
            [None, 0.5, 0.25],
        )

    def test_snapshot_loader_does_not_require_legacy_compatibility_file(self):
        published = self.publish_single(
            'Single_0000000001\t1\t23\tACGT/1/0\t1\t+\n'
        )
        snapshot = self.temp_dir / 'snapshot'
        snapshot.mkdir()
        copied = {}
        for key in (
            'normalized_genes',
            'normalized_guides',
            'normalized_gene_guide_manifest',
        ):
            destination = snapshot / Path(published[key]).name
            shutil.copy2(published[key], destination)
            copied[key] = destination

        loaded = normalization.load_normalized_gene_guide_data(
            copied['normalized_genes'],
            copied['normalized_guides'],
            copied['normalized_gene_guide_manifest'],
        )

        self.assertEqual(len(loaded['genes']), 1)
        self.assertEqual(len(loaded['guides']), 1)
        self.assertFalse(
            (snapshot / normalization.LEGACY_FINAL_FILENAME).exists()
        )

    def test_corrupted_current_artifact_fails_closed(self):
        published = self.publish_single(
            'Single_0000000001\t1\t23\tACGT/1/0\t1\t+\n'
        )
        with Path(published['normalized_guides']).open(
            'a', encoding='utf-8'
        ) as handle:
            handle.write('tampered\n')

        with self.assertRaisesRegex(
            normalization.GeneGuideNormalizationError,
            'size mismatch',
        ):
            normalization.resolve_current_normalized_gene_guide_paths(
                self.temp_dir
            )

    def test_pointer_publication_failure_preserves_previous_generation(self):
        first = self.publish_single(
            'Single_0000000001\t1\t23\tACGT/1/0\t1\t+\n'
        )
        pointer_path = self.temp_dir / normalization.CURRENT_POINTER_FILENAME
        previous_pointer = pointer_path.read_bytes()
        original_writer = normalization._write_json_atomic

        def fail_pointer(path, payload):
            if Path(path).name == normalization.CURRENT_POINTER_FILENAME:
                raise OSError('simulated pointer publication failure')
            return original_writer(path, payload)

        with mock.patch.object(
            normalization, '_write_json_atomic', side_effect=fail_pointer
        ):
            with self.assertRaisesRegex(OSError, 'simulated pointer'):
                normalization.publish_normalized_gene_guide_generation(
                    self.annotation_path,
                    self.temp_dir,
                    single_bed_path=self.temp_dir / 'single.bed',
                )

        self.assertEqual(pointer_path.read_bytes(), previous_pointer)
        current = normalization.resolve_current_normalized_gene_guide_paths(
            self.temp_dir
        )
        self.assertEqual(
            Path(current['normalized_genes']),
            Path(first['normalized_genes']),
        )

    def test_changed_annotation_source_invalidates_current_generation(self):
        self.publish_single(
            'Single_0000000001\t1\t23\tACGT/1/0\t1\t+\n'
        )
        original = self.annotation_path.read_text(encoding='utf-8')
        self.annotation_path.write_text(
            original.replace('GENE1', 'GENE2'), encoding='utf-8'
        )

        with self.assertRaisesRegex(
            normalization.GeneGuideNormalizationError,
            'source is stale: annotation_classified',
        ):
            normalization.resolve_current_normalized_gene_guide_paths(
                self.temp_dir
            )

    def test_semantically_invalid_staging_data_is_not_published(self):
        self.write_annotation([self.annotation_row()])

        def write_invalid_guides(path, connection):
            del connection
            Path(path).write_text('wrong\theader\n', encoding='utf-8')
            return 0

        with mock.patch.object(
            normalization,
            '_write_guides_from_index',
            side_effect=write_invalid_guides,
        ):
            with self.assertRaisesRegex(
                normalization.GeneGuideNormalizationError,
                'Normalized guide columns do not match',
            ):
                normalization.publish_normalized_gene_guide_generation(
                    self.annotation_path, self.temp_dir
                )

        self.assertFalse(
            (self.temp_dir / normalization.CURRENT_POINTER_FILENAME).exists()
        )

    def test_conflicting_gene_identity_is_rejected(self):
        self.write_annotation([
            self.annotation_row(
                query_id='query-a',
                subject_id='subject-a',
                classification='Multi',
                unique_id='Multi_0000000001',
            ),
            self.annotation_row(
                query_id='query-b',
                subject_id='subject-b',
                classification='Multi',
                unique_id='Multi_0000000001',
            ),
        ])

        with self.assertRaisesRegex(
            normalization.GeneGuideNormalizationError,
            'Conflicting subject_id',
        ):
            normalization.publish_normalized_gene_guide_generation(
                self.annotation_path, self.temp_dir
            )
        self.assertFalse(
            (self.temp_dir / normalization.CURRENT_POINTER_FILENAME).exists()
        )

    def test_resolver_returns_empty_mapping_for_legacy_runs(self):
        self.assertEqual(
            normalization.resolve_current_normalized_gene_guide_paths(
                self.temp_dir
            ),
            {},
        )

    def test_broken_current_pointer_symlink_fails_closed(self):
        pointer_path = self.temp_dir / normalization.CURRENT_POINTER_FILENAME
        pointer_path.symlink_to('missing-pointer.json')

        with self.assertRaisesRegex(
            normalization.GeneGuideNormalizationError,
            'Invalid normalized generation pointer',
        ):
            normalization.resolve_current_normalized_gene_guide_paths(
                self.temp_dir
            )

    def test_legacy_projection_lookup_uses_gene_ordinal_index(self):
        self.write_annotation([self.annotation_row()])
        _, rows = normalization._read_annotation(self.annotation_path)
        genes = normalization._canonical_genes(rows)
        bed_path = self.temp_dir / 'single.bed'
        bed_path.write_text(
            'Single_0000000001\t1\t23\tACGT/1/0\t1\t+\n',
            encoding='utf-8',
        )
        connection = normalization._create_guide_index(
            self.temp_dir / 'guide_index.sqlite3',
            None,
            bed_path,
            genes,
        )
        try:
            plan = connection.execute(
                'EXPLAIN QUERY PLAN SELECT raw_guide_target FROM guides '
                'WHERE unique_id = ? ORDER BY first_ordinal',
                ('Single_0000000001',),
            ).fetchall()
        finally:
            connection.close()

        self.assertIn(
            'guides_by_gene_ordinal',
            ' '.join(str(column) for row in plan for column in row),
        )


if __name__ == '__main__':
    unittest.main()
