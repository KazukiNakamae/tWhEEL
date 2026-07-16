import csv
import io
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

try:
    import requests  # noqa: F401
except ModuleNotFoundError:
    class RequestException(Exception):
        pass

    class ConnectionError(RequestException):
        pass

    class Timeout(RequestException):
        pass

    class HTTPError(RequestException):
        pass

    requests_module = types.ModuleType('requests')
    requests_module.exceptions = types.SimpleNamespace(
        RequestException=RequestException,
        ConnectionError=ConnectionError,
        Timeout=Timeout,
        HTTPError=HTTPError,
    )
    sys.modules['requests'] = requests_module

from src.docker.add_go_columns import add_go_columns


class AddGoColumnsTests(unittest.TestCase):
    @staticmethod
    def _empty_annotation():
        return {
            'GO_BP': '',
            'GO_MF': '',
            'GO_CC': '',
            'GO_ALL': '',
        }

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

    def test_transient_fetch_failure_is_retried(self):
        response = mock.Mock()
        response.text = (
            'Entry\tGene Ontology (biological process)\t'
            'Gene Ontology (molecular function)\t'
            'Gene Ontology (cellular component)\n'
            'P04637\tDNA damage response [GO:0006974]\t\t\n'
        )
        response.raise_for_status.return_value = None
        request_error = add_go_columns.requests.exceptions.RequestException(
            'temporary network failure'
        )

        with mock.patch.object(
            add_go_columns.requests,
            'get',
            side_effect=[request_error, response],
            create=True,
        ) as get_mock, mock.patch.object(
            add_go_columns.time, 'sleep'
        ) as sleep_mock:
            result = add_go_columns.fetch_go_batch(['P04637'])

        self.assertEqual(get_mock.call_count, 2)
        sleep_mock.assert_called_once_with(1)
        self.assertEqual(result['P04637']['GO_ALL'], 'GO:0006974')

    def test_repeated_fetch_failure_is_not_missing_go_annotation(self):
        request_error = add_go_columns.requests.exceptions.ConnectionError(
            'persistent network failure'
        )

        with mock.patch.object(
            add_go_columns.requests,
            'get',
            side_effect=request_error,
            create=True,
        ) as get_mock, mock.patch.object(
            add_go_columns.time, 'sleep'
        ) as sleep_mock:
            with self.assertRaises(
                add_go_columns.UniProtAPIError
            ) as error_context:
                add_go_columns.fetch_go_batch(['P04637'])

        self.assertEqual(get_mock.call_count, add_go_columns.API_MAX_ATTEMPTS)
        self.assertEqual(
            [call.args for call in sleep_mock.call_args_list],
            [(1,), (2,)],
        )
        error_message = str(error_context.exception)
        self.assertIn('persistent network failure', error_message)
        self.assertIn(
            'Network access is required for UniProt GO annotation',
            error_message,
        )
        self.assertIn('connect to the internet', error_message)
        self.assertIn('rerun the tWhEEL go_analysis step', error_message)

    def test_repeated_timeout_includes_network_and_rerun_guidance(self):
        request_error = add_go_columns.requests.exceptions.Timeout(
            'request timed out'
        )

        with mock.patch.object(
            add_go_columns.requests,
            'get',
            side_effect=request_error,
            create=True,
        ), mock.patch.object(add_go_columns.time, 'sleep'):
            with self.assertRaises(
                add_go_columns.UniProtAPIError
            ) as error_context:
                add_go_columns.fetch_go_batch(['P04637'])

        error_message = str(error_context.exception)
        self.assertIn(
            'Network access is required for UniProt GO annotation',
            error_message,
        )
        self.assertIn('connect to the internet', error_message)
        self.assertIn('rerun the tWhEEL go_analysis step', error_message)

    def test_repeated_http_error_does_not_claim_network_is_unavailable(self):
        request_error = add_go_columns.requests.exceptions.HTTPError(
            '503 Service Unavailable'
        )

        with mock.patch.object(
            add_go_columns.requests,
            'get',
            side_effect=request_error,
            create=True,
        ), mock.patch.object(add_go_columns.time, 'sleep'):
            with self.assertRaises(
                add_go_columns.UniProtAPIError
            ) as error_context:
                add_go_columns.fetch_go_batch(['P04637'])

        error_message = str(error_context.exception)
        self.assertIn('503 Service Unavailable', error_message)
        self.assertNotIn(
            'Network access is required for UniProt GO annotation',
            error_message,
        )
        self.assertNotIn('connect to the internet', error_message)

    def test_cli_reports_network_guidance_to_stderr_and_exits_one(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'input.tsv'
            output_path = temp_dir / 'output.tsv'
            input_path.write_text(
                'query_id\tsubject_id\nquery_1\tP04637\n',
                encoding='utf-8',
            )
            request_error = add_go_columns.requests.exceptions.Timeout(
                'request timed out'
            )
            stderr = io.StringIO()

            with mock.patch.object(
                add_go_columns.requests,
                'get',
                side_effect=request_error,
                create=True,
            ), mock.patch.object(
                add_go_columns.time, 'sleep'
            ), mock.patch.object(
                add_go_columns.sys,
                'argv',
                ['add_go_columns.py', str(input_path), str(output_path)],
            ), mock.patch.object(add_go_columns.sys, 'stderr', stderr):
                with self.assertRaises(SystemExit) as error_context:
                    add_go_columns.main()

            self.assertEqual(error_context.exception.code, 1)
            error_output = stderr.getvalue()
            self.assertIn(
                'Network access is required for UniProt GO annotation',
                error_output,
            )
            self.assertIn('connect to the internet', error_output)
            self.assertIn(
                'rerun the tWhEEL go_analysis step',
                error_output,
            )

    def test_malformed_success_response_is_not_missing_go_annotation(self):
        response = mock.Mock()
        response.text = '<html>unexpected proxy response</html>'
        response.raise_for_status.return_value = None

        with mock.patch.object(
            add_go_columns.requests,
            'get',
            return_value=response,
            create=True,
        ):
            with self.assertRaisesRegex(
                add_go_columns.UniProtAPIError,
                'invalid TSV response',
            ):
                add_go_columns.fetch_go_batch(['P04637'])

    def test_partial_batch_failure_preserves_previous_complete_output(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'input.tsv'
            output_path = temp_dir / 'output.tsv'

            with input_path.open('w', encoding='utf-8', newline='') as handle:
                writer = csv.writer(handle, delimiter='\t')
                writer.writerow(['query_id', 'subject_id'])
                for index in range(add_go_columns.BATCH_SIZE + 1):
                    writer.writerow([f'query_{index}', f'P{index:05d}'])

            previous_output = 'previous complete output\n'
            output_path.write_text(previous_output, encoding='utf-8')

            call_count = 0

            def fetch_with_second_batch_failure(uniprot_ids):
                nonlocal call_count
                call_count += 1
                if call_count == 2:
                    raise add_go_columns.UniProtAPIError(
                        'second batch failed'
                    )
                return {
                    accession: self._empty_annotation()
                    for accession in uniprot_ids
                }

            with mock.patch.object(
                add_go_columns,
                'fetch_go_batch',
                side_effect=fetch_with_second_batch_failure,
            ), mock.patch.object(
                add_go_columns.time, 'sleep'
            ), mock.patch.object(
                add_go_columns.sys,
                'argv',
                ['add_go_columns.py', str(input_path), str(output_path)],
            ):
                with self.assertRaises(SystemExit) as error_context:
                    add_go_columns.main()

            self.assertEqual(call_count, 2)
            self.assertEqual(error_context.exception.code, 1)
            self.assertEqual(
                output_path.read_text(encoding='utf-8'),
                previous_output,
            )
            self.assertEqual(
                list(temp_dir.glob(f'.{output_path.name}.*.tmp')),
                [],
            )

    def test_success_atomically_replaces_output(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'input.tsv'
            output_path = temp_dir / 'output.tsv'
            input_path.write_text(
                'query_id\tsubject_id\nquery_1\tsp|P04637|P53_HUMAN\n',
                encoding='utf-8',
            )
            output_path.write_text('old output\n', encoding='utf-8')
            annotation = {
                'GO_BP': 'GO:0006974',
                'GO_MF': '',
                'GO_CC': '',
                'GO_ALL': 'GO:0006974',
            }

            with mock.patch.object(
                add_go_columns,
                'fetch_go_batch',
                return_value={'P04637': annotation},
            ), mock.patch.object(add_go_columns.time, 'sleep'):
                processed_rows = add_go_columns.add_go_columns(
                    input_path, output_path
                )

            self.assertEqual(processed_rows, 1)
            with output_path.open(encoding='utf-8', newline='') as handle:
                rows = list(csv.DictReader(handle, delimiter='\t'))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]['GO_ALL'], 'GO:0006974')

    def test_subject_id_is_resolved_by_header_name_for_canonical_gene_rows(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'canonical_genes.tsv'
            output_path = temp_dir / 'canonical_genes_go.tsv'
            input_path.write_text(
                'unique_id\tguide_count\tclassification\tsubject_id\n'
                'Single_1\t2\tSingle\tsp|P04637|P53_HUMAN\n',
                encoding='utf-8',
            )
            annotation = {
                'GO_BP': 'GO:0006974',
                'GO_MF': '',
                'GO_CC': '',
                'GO_ALL': 'GO:0006974',
            }

            with mock.patch.object(
                add_go_columns,
                'fetch_go_batch',
                return_value={'P04637': annotation},
            ) as fetch_mock, mock.patch.object(add_go_columns.time, 'sleep'):
                add_go_columns.add_go_columns(input_path, output_path)

            fetch_mock.assert_called_once_with(['P04637'])
            with output_path.open(encoding='utf-8', newline='') as handle:
                reader = csv.DictReader(handle, delimiter='\t')
                rows = list(reader)
            self.assertEqual(reader.fieldnames, [
                'unique_id', 'guide_count', 'classification', 'subject_id',
                'GO_BP', 'GO_MF', 'GO_CC', 'GO_ALL',
            ])
            self.assertEqual(rows[0]['unique_id'], 'Single_1')
            self.assertEqual(rows[0]['guide_count'], '2')
            self.assertEqual(rows[0]['GO_ALL'], 'GO:0006974')

    def test_missing_subject_id_header_preserves_previous_output(self):
        with tempfile.TemporaryDirectory() as temp_dir_name:
            temp_dir = Path(temp_dir_name)
            input_path = temp_dir / 'input.tsv'
            output_path = temp_dir / 'output.tsv'
            input_path.write_text(
                'unique_id\tguide_count\nSingle_1\t2\n',
                encoding='utf-8',
            )
            output_path.write_text('previous output\n', encoding='utf-8')

            with mock.patch.object(
                add_go_columns, 'fetch_go_batch'
            ) as fetch_mock:
                with self.assertRaisesRegex(ValueError, "'subject_id' column"):
                    add_go_columns.add_go_columns(input_path, output_path)

            fetch_mock.assert_not_called()
            self.assertEqual(
                output_path.read_text(encoding='utf-8'), 'previous output\n'
            )
            self.assertEqual(
                list(temp_dir.glob(f'.{output_path.name}.*.tmp')), []
            )


if __name__ == '__main__':
    unittest.main()
