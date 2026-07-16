import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from src.modules.annotation import Annotation
from src.modules.sgd_annotation import Sgd_annotation


class MultiGeneAlignmentContractTests(unittest.TestCase):
    def test_annotation_rejects_success_without_expected_output(self):
        self._assert_empty_output_is_rejected(
            Annotation,
            'src.modules.annotation.subprocess.run',
        )

    def test_sgd_annotation_rejects_success_without_expected_output(self):
        self._assert_empty_output_is_rejected(
            Sgd_annotation,
            'src.modules.sgd_annotation.subprocess.run',
        )

    def _assert_empty_output_is_rejected(self, module_class, run_target):
        with tempfile.TemporaryDirectory() as temp_dir:
            multi_genes_dir = Path(temp_dir) / 'multi_genes'
            multi_genes_dir.mkdir()
            (multi_genes_dir / 'Multi_1.fa').write_text(
                '>Multi_1_tx1\nACGT\n', encoding='utf-8'
            )
            aligned_dir = Path(temp_dir) / 'multi_genes_aligned'
            aligned_dir.mkdir()
            stale_output = aligned_dir / 'Multi_1_aligned_consensus.fa'
            stale_output.write_text(
                '>Multi_1\nSTALE\n', encoding='utf-8'
            )
            module = module_class({
                'output': {'intermediate': temp_dir},
            })
            successful_result = subprocess.CompletedProcess(
                args=[], returncode=0, stdout='', stderr=''
            )

            with patch(run_target, return_value=successful_result):
                with self.assertRaisesRegex(
                    RuntimeError, 'missing=Multi_1'
                ):
                    module._run_alignment_process({'Multi_1'})
            self.assertFalse(stale_output.exists())


if __name__ == '__main__':
    unittest.main()
