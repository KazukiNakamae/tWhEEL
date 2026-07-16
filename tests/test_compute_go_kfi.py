import unittest
from unittest import mock

import pandas as pd

from src.docker.compute_go_kfi.compute_go_kfi import (
    calculate_go_kfi,
    validate_input_data,
)


class CalculateGoKfiTests(unittest.TestCase):
    def test_duplicate_annotation_rows_count_each_gene_once(self):
        data = pd.DataFrame([
            {
                'unique_id': 'Multi_1',
                'guide_targets': 'guide_a,guide_b',
                'GO_ALL': 'GO:0000001;GO:0000001',
            },
            {
                'unique_id': 'Multi_1',
                'guide_targets': 'guide_a,guide_b',
                'GO_ALL': 'GO:0000001',
            },
            {
                'unique_id': 'Single_1',
                'guide_targets': '',
                'GO_ALL': 'GO:0000001',
            },
        ])

        output = calculate_go_kfi(data)
        self.assertEqual(list(output.columns), [
            'GO_ID',
            'total_genes',
            'genes_with_guides',
            'total_guides',
            'kfi_score',
            'targetable_type',
        ])
        result = output.set_index('GO_ID').loc['GO:0000001']

        self.assertEqual(result['total_genes'], 2)
        self.assertEqual(result['genes_with_guides'], 1)
        self.assertEqual(result['total_guides'], 2)
        self.assertEqual(result['kfi_score'], 0.5)
        self.assertEqual(result['targetable_type'], 'Partially')

    def test_gene_uses_union_of_go_terms_and_guides_across_rows(self):
        data = pd.DataFrame([
            {
                'unique_id': 'Multi_1',
                'guide_targets': 'guide_a',
                'GO_ALL': 'GO:0000001',
            },
            {
                'unique_id': 'Multi_1',
                'guide_targets': 'guide_b',
                'GO_ALL': 'GO:0000002',
            },
            {
                'unique_id': 'Single_1',
                'guide_targets': '',
                'GO_ALL': 'GO:0000001',
            },
        ])

        result = calculate_go_kfi(data).set_index('GO_ID')

        self.assertEqual(result.loc['GO:0000001', 'total_genes'], 2)
        self.assertEqual(result.loc['GO:0000001', 'genes_with_guides'], 1)
        self.assertEqual(result.loc['GO:0000001', 'total_guides'], 2)
        self.assertEqual(result.loc['GO:0000001', 'kfi_score'], 0.5)
        self.assertEqual(result.loc['GO:0000002', 'total_genes'], 1)
        self.assertEqual(result.loc['GO:0000002', 'genes_with_guides'], 1)
        self.assertEqual(result.loc['GO:0000002', 'total_guides'], 2)
        self.assertEqual(result.loc['GO:0000002', 'kfi_score'], 1.0)
        self.assertEqual(result.loc['GO:0000002', 'targetable_type'], 'Fully')

    def test_missing_unique_id_is_rejected(self):
        data = pd.DataFrame([
            {
                'unique_id': pd.NA,
                'guide_targets': 'guide_a',
                'GO_ALL': 'GO:0000001',
            },
        ])

        with self.assertRaisesRegex(ValueError, 'valid unique_id'):
            calculate_go_kfi(data)

    def test_canonical_guide_count_is_preferred_without_parsing_targets(self):
        data = pd.DataFrame([
            {
                'unique_id': 'Single_1',
                'subject_id': 'sp|P04637|P53_HUMAN',
                'guide_count': 3,
                'guide_targets': 'not/a/serialized/legacy/field',
                'GO_ALL': 'GO:0000001;GO:0000002',
            },
            {
                'unique_id': 'Single_2',
                'subject_id': 'sp|P38398|BRCA1_HUMAN',
                'guide_count': 0,
                'guide_targets': 'legacy_value_that_must_not_be_counted',
                'GO_ALL': 'GO:0000001',
            },
        ])

        with mock.patch(
            'src.docker.compute_go_kfi.compute_go_kfi.'
            '_parse_guide_target_set',
            side_effect=AssertionError('legacy parser was called'),
        ):
            result = calculate_go_kfi(data).set_index('GO_ID')

        self.assertEqual(result.loc['GO:0000001', 'total_genes'], 2)
        self.assertEqual(result.loc['GO:0000001', 'genes_with_guides'], 1)
        self.assertEqual(result.loc['GO:0000001', 'total_guides'], 3)
        self.assertEqual(result.loc['GO:0000001', 'kfi_score'], 0.5)
        self.assertEqual(result.loc['GO:0000002', 'total_genes'], 1)
        self.assertEqual(result.loc['GO:0000002', 'total_guides'], 3)
        self.assertEqual(result.loc['GO:0000002', 'targetable_type'], 'Fully')
        self.assertTrue(validate_input_data(data))

    def test_conflicting_guide_counts_for_one_gene_are_rejected(self):
        data = pd.DataFrame([
            {
                'unique_id': 'Multi_1',
                'guide_count': 1,
                'GO_ALL': 'GO:0000001',
            },
            {
                'unique_id': 'Multi_1',
                'guide_count': 2,
                'GO_ALL': 'GO:0000002',
            },
        ])

        with self.assertRaisesRegex(
            ValueError, "Conflicting guide_count values for 'Multi_1'"
        ):
            calculate_go_kfi(data)
        self.assertFalse(validate_input_data(data))

    def test_invalid_canonical_guide_counts_are_rejected(self):
        invalid_counts = [-1, 1.5, 'two', '', pd.NA]
        for guide_count in invalid_counts:
            with self.subTest(guide_count=guide_count):
                data = pd.DataFrame([{
                    'unique_id': 'Single_1',
                    'guide_count': guide_count,
                    'GO_ALL': 'GO:0000001',
                }])

                with self.assertRaisesRegex(
                    ValueError, 'must be a non-negative integer'
                ):
                    calculate_go_kfi(data)

    def test_repeated_identical_guide_count_preserves_go_union_semantics(self):
        data = pd.DataFrame([
            {
                'unique_id': 'Multi_1',
                'guide_count': 2,
                'GO_ALL': 'GO:0000001',
            },
            {
                'unique_id': 'Multi_1',
                'guide_count': 2,
                'GO_ALL': 'GO:0000002',
            },
        ])

        result = calculate_go_kfi(data).set_index('GO_ID')

        for go_id in ('GO:0000001', 'GO:0000002'):
            self.assertEqual(result.loc[go_id, 'total_genes'], 1)
            self.assertEqual(result.loc[go_id, 'genes_with_guides'], 1)
            self.assertEqual(result.loc[go_id, 'total_guides'], 2)
            self.assertEqual(result.loc[go_id, 'targetable_type'], 'Fully')


if __name__ == '__main__':
    unittest.main()
