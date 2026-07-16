import csv
import tempfile
import unittest
from pathlib import Path

from src.modules.annotation import Annotation


class AnnotationParserTests(unittest.TestCase):
    def test_classifies_custom_subject_without_description(self):
        ggsearch_output = """\
1>>>query_only - 120 aa
The best scores are: n-w bits E(1)
custom_protein_42 ( 120) 795 105.4 0

>>>query_only, 120 aa vs protein_db.fasta library
>>custom_protein_42
>>><<<
"""

        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'ggsearch_results.txt'
            output_path = temp_dir / 'annotation_table.txt'
            input_path.write_text(ggsearch_output, encoding='utf-8')
            annotation = Annotation({
                'output': {'intermediate': str(temp_dir)},
            })

            annotation._parse_ggsearch_output(input_path, output_path)
            classified_path = Path(annotation.classify_genes(output_path))

            with classified_path.open(
                encoding='utf-8', newline=''
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['query_id'], 'query_only')
        self.assertEqual(rows[0]['subject_id'], 'custom_protein_42')
        self.assertEqual(rows[0]['description'], '')
        self.assertEqual(rows[0]['gene_name'], 'unknown')
        self.assertEqual(rows[0]['classification'], 'Single')
        self.assertEqual(rows[0]['unique_id'], 'Single_0000000001')

    def test_accepts_query_header_without_fasta_description(self):
        ggsearch_output = """\
1>>>query_only - 120 aa
The best scores are: n-w bits E(1)
sp|P00001|TEST Example protein GN=TEST ( 120) 795 105.4 0

>>>query_only, 120 aa vs protein_db.fasta library
>>sp|P00001|TEST Example protein GN=TEST
>>><<<
"""

        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'ggsearch_results.txt'
            output_path = temp_dir / 'annotation_table.txt'
            input_path.write_text(ggsearch_output, encoding='utf-8')
            annotation = Annotation({
                'output': {'intermediate': str(temp_dir)},
            })

            annotation._parse_ggsearch_output(input_path, output_path)

            with output_path.open(encoding='utf-8', newline='') as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]['query_id'], 'query_only')
        self.assertEqual(rows[0]['subject_id'], 'sp|P00001|TEST')
        self.assertEqual(rows[0]['description'], 'Example protein GN=TEST')
        self.assertEqual(rows[0]['length'], '120')
        self.assertEqual(rows[0]['bit_score'], '105.4')
        self.assertEqual(rows[0]['evalue'], '0')

    def test_accepts_swissprot_trembl_ncbi_and_custom_subject_ids(self):
        cases = [
            {
                'query_id': 'query_sp',
                'subject_id': 'sp|P17019|ZN708_HUMAN',
                'summary': 'Zinc finger protein 708 OS=H',
                'full': 'Zinc finger protein 708 OS=Homo sapiens GN=ZNF708',
            },
            {
                'query_id': 'query_tr',
                'subject_id': 'tr|A0A123|A0A123_MOUSE',
                'summary': 'Unreviewed protein OS=M',
                'full': 'Unreviewed protein OS=Mus musculus GN=TrGene',
            },
            {
                'query_id': 'query_ncbi',
                'subject_id': 'ref|XP_012345.1|',
                'summary': 'Predicted protein OS=H',
                'full': 'Predicted protein OS=Homo sapiens GN=NCBIGene',
                'raw_score': '-20',
            },
            {
                'query_id': 'query_custom',
                'subject_id': 'custom_protein_42',
                'summary': 'Locally curated protein',
                'full': 'Locally curated protein GN=CustomGene source=internal',
            },
        ]

        standard_sections = []
        m10_sections = []
        for index, case in enumerate(cases, 1):
            standard_sections.extend([
                f"{index}>>>{case['query_id']} test protein - 100 aa",
                'The best scores are: n-w bits E(1)',
                f"{case['subject_id']} {case['summary']} ( 100) "
                f"{case.get('raw_score', '200')} 50.0 1e-20",
                'sp|WRONG|SECOND_HIT This must not replace the best hit ( 90) 180 45.0 1e-10',
                '',
            ])
            m10_sections.extend([
                f">>>{case['query_id']}, 100 aa vs protein_db.fasta library",
                f">>{case['subject_id']} {case['full']}",
                '>>><<<',
            ])

        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'ggsearch_results.txt'
            output_path = temp_dir / 'annotation_table.txt'
            input_path.write_text(
                '\n'.join(standard_sections + m10_sections),
                encoding='utf-8',
            )
            annotation = Annotation({
                'output': {'intermediate': str(temp_dir)},
            })

            annotation._parse_ggsearch_output(input_path, output_path)

            with output_path.open(encoding='utf-8', newline='') as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))

        self.assertEqual(len(rows), len(cases))
        self.assertEqual(
            [row['subject_id'] for row in rows],
            [case['subject_id'] for case in cases],
        )
        self.assertEqual(
            [row['description'] for row in rows],
            [case['full'] for case in cases],
        )
        self.assertTrue(all(row['length'] == '100' for row in rows))
        self.assertTrue(all(row['bit_score'] == '50.0' for row in rows))
        self.assertTrue(all(row['evalue'] == '1e-20' for row in rows))


if __name__ == '__main__':
    unittest.main()
