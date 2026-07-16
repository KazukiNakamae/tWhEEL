import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.docker.run_alignment_and_mask.run_alignment_and_mask import (
    AlignmentProcessor,
)


class AlignmentProcessorTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        temp_path = Path(self.temp_dir.name)
        self.processor = AlignmentProcessor(
            temp_path / 'input',
            temp_path / 'output',
        )
        self.sequences = [
            SeqRecord(Seq('ACGT'), id='Multi_1_transcript_1'),
            SeqRecord(Seq('ACT'), id='Multi_1_transcript_2'),
        ]

    @mock.patch(
        'src.docker.run_alignment_and_mask.run_alignment_and_mask.subprocess.run'
    )
    def test_mafft_nonzero_exit_raises(self, run_mock):
        run_mock.return_value = subprocess.CompletedProcess(
            args=['mafft'],
            returncode=42,
            stderr='simulated MAFFT failure',
        )

        with self.assertRaisesRegex(RuntimeError, 'exit code 42'):
            self.processor.run_mafft_alignment(self.sequences, 'Multi_1')

    @mock.patch(
        'src.docker.run_alignment_and_mask.run_alignment_and_mask.subprocess.run'
    )
    def test_empty_mafft_output_raises(self, run_mock):
        run_mock.return_value = subprocess.CompletedProcess(
            args=['mafft'],
            returncode=0,
            stderr='',
        )

        with self.assertRaisesRegex(RuntimeError, 'produced no aligned sequences'):
            self.processor.run_mafft_alignment(self.sequences, 'Multi_1')

    @mock.patch(
        'src.docker.run_alignment_and_mask.run_alignment_and_mask.subprocess.run'
    )
    def test_incomplete_mafft_output_raises(self, run_mock):
        def write_incomplete_alignment(*args, **kwargs):
            kwargs['stdout'].write('>Multi_1_transcript_1\nACGT\n')
            return subprocess.CompletedProcess(
                args=args[0],
                returncode=0,
                stderr='',
            )

        run_mock.side_effect = write_incomplete_alignment

        with self.assertRaisesRegex(RuntimeError, 'returned 1 of 2 sequences'):
            self.processor.run_mafft_alignment(self.sequences, 'Multi_1')

    @mock.patch(
        'src.docker.run_alignment_and_mask.run_alignment_and_mask.subprocess.run'
    )
    def test_valid_mafft_output_can_generate_consensus(self, run_mock):
        def write_alignment(*args, **kwargs):
            kwargs['stdout'].write(
                '>Multi_1_transcript_1\nAC-G\n'
                '>Multi_1_transcript_2\nACTG\n'
            )
            return subprocess.CompletedProcess(
                args=args[0],
                returncode=0,
                stderr='',
            )

        run_mock.side_effect = write_alignment

        aligned = self.processor.run_mafft_alignment(
            self.sequences,
            'Multi_1',
        )
        consensus = self.processor.generate_consensus(aligned, 'Multi_1')

        self.assertEqual(len(aligned), 2)
        self.assertEqual({len(record.seq) for record in aligned}, {4})
        self.assertEqual(str(consensus.seq), 'ACTG')

    def test_consensus_rejects_unaligned_sequences(self):
        with self.assertRaisesRegex(ValueError, 'not a valid alignment'):
            self.processor.generate_consensus(self.sequences, 'Multi_1')

    def test_run_raises_after_any_group_failure(self):
        groups = {
            'Multi_1': self.sequences,
            'Multi_2': self.sequences,
        }
        self.processor.group_sequences_by_multi_id = mock.Mock(
            return_value=groups,
        )
        self.processor.process_multi_gene_group = mock.Mock(
            side_effect=[Path('Multi_1.fa'), RuntimeError('MAFFT failed')],
        )

        with self.assertRaisesRegex(RuntimeError, '1 of 2 groups'):
            self.processor.run()

        self.assertEqual(
            self.processor.process_multi_gene_group.call_count,
            2,
        )


if __name__ == '__main__':
    unittest.main()
