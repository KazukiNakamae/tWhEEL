import unittest

import pandas as pd

from src.docker.compute_go_kfi.compute_go_kfi import calculate_go_kfi


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


if __name__ == '__main__':
    unittest.main()
