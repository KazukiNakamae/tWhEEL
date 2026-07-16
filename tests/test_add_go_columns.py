import sys
import types
import unittest
from unittest import mock

try:
    import requests  # noqa: F401
except ModuleNotFoundError:
    requests_module = types.ModuleType('requests')
    requests_module.exceptions = types.SimpleNamespace(RequestException=Exception)
    sys.modules['requests'] = requests_module

from src.docker.add_go_columns import add_go_columns


class AddGoColumnsTests(unittest.TestCase):
    def test_classifies_terms_using_uniprot_aspects(self):
        classified = add_go_columns.classify_go_terms({
            'GO_BP': ['GO:0006974'],
            'GO_MF': ['GO:0003677', 'GO:0003677'],
            'GO_CC': ['GO:0005634'],
        })

        self.assertEqual(classified['GO_BP'], 'GO:0006974')
        self.assertEqual(classified['GO_MF'], 'GO:0003677')
        self.assertEqual(classified['GO_CC'], 'GO:0005634')
        self.assertEqual(
            classified['GO_ALL'],
            'GO:0006974;GO:0003677;GO:0005634',
        )
        self.assertNotIn('GO:0003677', classified['GO_BP'])

    def test_fetch_requests_and_parses_aspect_specific_fields(self):
        response = mock.Mock()
        response.text = (
            'Entry\tGene Ontology (biological process)\t'
            'Gene Ontology (molecular function)\t'
            'Gene Ontology (cellular component)\n'
            'P04637\tDNA damage response [GO:0006974]\t'
            'DNA binding [GO:0003677]; DNA binding [GO:0003677]\t'
            'nucleus [GO:0005634]\n'
        )
        response.raise_for_status.return_value = None

        with mock.patch.object(
            add_go_columns.requests,
            'get',
            return_value=response,
            create=True,
        ) as get_mock:
            result = add_go_columns.fetch_go_batch(['P04637', 'NO_GO_ENTRY'])

        params = get_mock.call_args.kwargs['params']
        self.assertEqual(params['fields'], 'accession,go_p,go_f,go_c')
        self.assertEqual(result['P04637'], {
            'GO_BP': 'GO:0006974',
            'GO_MF': 'GO:0003677',
            'GO_CC': 'GO:0005634',
            'GO_ALL': 'GO:0006974;GO:0003677;GO:0005634',
        })
        self.assertEqual(result['NO_GO_ENTRY'], {
            'GO_BP': '',
            'GO_MF': '',
            'GO_CC': '',
            'GO_ALL': '',
        })


if __name__ == '__main__':
    unittest.main()
