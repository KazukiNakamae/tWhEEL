import csv
import tempfile
import unittest
from pathlib import Path

from src.modules.sgd_annotation import Sgd_annotation


class SgdAnnotationParserTests(unittest.TestCase):
    def test_preserves_gene_symbol_and_score_columns(self):
        ggsearch_output = """\
1>>>EIW12331 test protein - 383 aa
The best scores are: n-w bits E(1)
YAL060W|BDH1|SGDID:S000000056 NAD-dependent (R,R)- ( 383) 2584 154.8 0

"""

        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'ggsearch_results_sgd.txt'
            annotation_path = temp_dir / 'annotation_table.txt'
            input_path.write_text(ggsearch_output, encoding='utf-8')
            annotation = Sgd_annotation({
                'output': {'intermediate': str(temp_dir)},
            })

            annotation._parse_ggsearch_output(input_path, annotation_path)
            with annotation_path.open(encoding='utf-8', newline='') as handle:
                parsed_rows = list(csv.DictReader(handle, delimiter='\t'))

            classified_path = Path(annotation.classify_genes(annotation_path))
            with classified_path.open(encoding='utf-8', newline='') as handle:
                classified_rows = list(csv.DictReader(handle, delimiter='\t'))

        self.assertEqual(len(parsed_rows), 1)
        self.assertEqual(
            parsed_rows[0],
            {
                'query_id': 'EIW12331',
                'subject_id': 'SGDID:S000000056',
                'identity': '',
                'length': '383',
                'evalue': '0',
                'bit_score': '154.8',
                'description': 'NAD-dependent (R,R)-',
                'gene_name': 'BDH1',
                'annotation_source_database': 'SGD',
                'sgd_id': 'S000000056',
                'ensembl_id': 'YAL060W',
            },
        )
        self.assertEqual(len(classified_rows), 1)
        self.assertEqual(classified_rows[0]['gene_name'], 'BDH1')
        self.assertNotEqual(classified_rows[0]['gene_name'], 'NAD-dependent')
        self.assertEqual(classified_rows[0]['length'], '383')
        self.assertEqual(classified_rows[0]['bit_score'], '154.8')
        self.assertEqual(classified_rows[0]['evalue'], '0')
        self.assertEqual(classified_rows[0]['annotation_source_database'], 'SGD')
        self.assertEqual(classified_rows[0]['sgd_id'], 'S000000056')
        self.assertEqual(classified_rows[0]['ensembl_id'], 'YAL060W')
        self.assertEqual(len(classified_rows[0]), 13)


if __name__ == '__main__':
    unittest.main()
