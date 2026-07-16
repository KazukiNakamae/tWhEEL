import json
import os
import subprocess
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.docker.run_alignment_and_mask.run_alignment_and_mask import (
    AlignmentProcessor,
    GroupDescriptor,
)
from src.modules.multi_gene_alignment import annotation_workspace_lock


class AlignmentProcessorTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.temp_path = Path(self.temp_dir.name)
        self.input_dir = self.temp_path / "input"
        self.output_dir = self.temp_path / "output"
        self.cache_dir = self.temp_path / "cache"
        self.input_dir.mkdir()
        self.processor = AlignmentProcessor(
            self.input_dir,
            self.output_dir,
            cache_dir=self.cache_dir,
        )
        self.sequences = [
            SeqRecord(Seq("ACGT"), id="Multi_1_transcript_1"),
            SeqRecord(Seq("ACT"), id="Multi_1_transcript_2"),
        ]

    @staticmethod
    def _consensus(multi_id, sequence="ACTG"):
        return SeqRecord(
            Seq(sequence),
            id=multi_id,
            description="test consensus",
        )

    def _set_groups(self, groups):
        self.processor.group_sequences_by_multi_id = mock.Mock(
            return_value=groups
        )

    @mock.patch(
        "src.docker.run_alignment_and_mask.run_alignment_and_mask.subprocess.run"
    )
    def test_mafft_nonzero_exit_raises(self, run_mock):
        run_mock.return_value = subprocess.CompletedProcess(
            args=["mafft"],
            returncode=42,
            stderr="simulated MAFFT failure",
        )

        with self.assertRaisesRegex(RuntimeError, "exit code 42"):
            self.processor.run_mafft_alignment(self.sequences, "Multi_1")

    @mock.patch(
        "src.docker.run_alignment_and_mask.run_alignment_and_mask.subprocess.run"
    )
    def test_empty_mafft_output_raises(self, run_mock):
        run_mock.return_value = subprocess.CompletedProcess(
            args=["mafft"],
            returncode=0,
            stderr="",
        )

        with self.assertRaisesRegex(RuntimeError, "produced no aligned sequences"):
            self.processor.run_mafft_alignment(self.sequences, "Multi_1")

    @mock.patch(
        "src.docker.run_alignment_and_mask.run_alignment_and_mask.subprocess.run"
    )
    def test_incomplete_mafft_output_raises(self, run_mock):
        def write_incomplete_alignment(*args, **kwargs):
            kwargs["stdout"].write(">Multi_1_transcript_1\nACGT\n")
            return subprocess.CompletedProcess(
                args=args[0],
                returncode=0,
                stderr="",
            )

        run_mock.side_effect = write_incomplete_alignment

        with self.assertRaisesRegex(RuntimeError, "returned 1 of 2 sequences"):
            self.processor.run_mafft_alignment(self.sequences, "Multi_1")

    @mock.patch(
        "src.docker.run_alignment_and_mask.run_alignment_and_mask.subprocess.run"
    )
    def test_valid_mafft_output_can_generate_consensus(self, run_mock):
        def write_alignment(*args, **kwargs):
            kwargs["stdout"].write(
                ">Multi_1_transcript_1\nAC-G\n"
                ">Multi_1_transcript_2\nACTG\n"
            )
            return subprocess.CompletedProcess(
                args=args[0],
                returncode=0,
                stderr="",
            )

        run_mock.side_effect = write_alignment

        aligned = self.processor.run_mafft_alignment(
            self.sequences,
            "Multi_1",
        )
        consensus = self.processor.generate_consensus(aligned, "Multi_1")

        command = run_mock.call_args.args[0]
        self.assertEqual(command[command.index("--thread") + 1], "1")
        self.assertEqual(len(aligned), 2)
        self.assertEqual({len(record.seq) for record in aligned}, {4})
        self.assertEqual(str(consensus.seq), "ACTG")

    def test_consensus_rejects_unaligned_sequences(self):
        with self.assertRaisesRegex(ValueError, "not a valid alignment"):
            self.processor.generate_consensus(self.sequences, "Multi_1")

    def test_standard_group_discovery_defers_parsing_to_worker(self):
        (self.input_dir / "Multi_1.fa").write_text(
            ">Multi_1_transcript_1\nACGT\n"
            ">Multi_1_transcript_2\nACT\n",
            encoding="utf-8",
        )
        (self.input_dir / "Multi_2.fa").write_text(
            ">Multi_2_transcript_1\nGGCC\n",
            encoding="utf-8",
        )

        with mock.patch(
            "src.docker.run_alignment_and_mask.run_alignment_and_mask.SeqIO.parse",
            wraps=SeqIO.parse,
        ) as parse_mock:
            groups = self.processor.group_sequences_by_multi_id()
            self.assertEqual(parse_mock.call_count, 0)
            loaded = self.processor._load_group_sequences(
                "Multi_1", groups["Multi_1"]
            )
            self.assertEqual(parse_mock.call_count, 1)

        self.assertEqual(set(groups), {"Multi_1", "Multi_2"})
        self.assertTrue(
            all(isinstance(group, GroupDescriptor) for group in groups.values())
        )
        self.assertTrue(
            all(group.strict_single_file for group in groups.values())
        )
        self.assertEqual(
            [record.id for record in loaded],
            ["Multi_1_transcript_1", "Multi_1_transcript_2"],
        )

    def test_nonstandard_group_file_uses_compatible_discovery_scan(self):
        (self.input_dir / "legacy_group.fa").write_text(
            ">Multi_1_transcript_1\nACGT\n"
            ">Multi_2_transcript_1\nGGCC\n",
            encoding="utf-8",
        )

        with mock.patch(
            "src.docker.run_alignment_and_mask.run_alignment_and_mask.SeqIO.parse",
            wraps=SeqIO.parse,
        ) as parse_mock:
            groups = self.processor.group_sequences_by_multi_id()
            first_group = self.processor._load_group_sequences(
                "Multi_1", groups["Multi_1"]
            )
            second_group = self.processor._load_group_sequences(
                "Multi_2", groups["Multi_2"]
            )

        self.assertEqual(set(groups), {"Multi_1", "Multi_2"})
        self.assertFalse(groups["Multi_1"].strict_single_file)
        self.assertEqual(parse_mock.call_count, 1)
        self.assertEqual(
            [record.id for record in first_group],
            ["Multi_1_transcript_1"],
        )
        self.assertEqual(
            [record.id for record in second_group],
            ["Multi_2_transcript_1"],
        )

    def test_default_cache_identity_is_derived_at_runtime(self):
        with mock.patch.object(
            AlignmentProcessor,
            "_automatic_tool_identity",
            return_value="auto-sha256:test-identity",
        ) as identity_mock:
            processor = AlignmentProcessor(
                self.input_dir,
                self.temp_path / "automatic-output",
                cache_dir=self.temp_path / "automatic-cache",
            )

        identity_mock.assert_called_once_with()
        self.assertEqual(
            processor.tool_identity,
            "auto-sha256:test-identity",
        )

    def test_skipped_lock_requires_and_verifies_parent_ownership(self):
        with self.assertRaisesRegex(ValueError, "parent_lock_token"):
            AlignmentProcessor(
                self.input_dir,
                self.output_dir,
                cache_enabled=False,
                acquire_run_lock=False,
            )

        lock_dir = self.temp_path / ".annotation-workspace.lock.d"
        lock_dir.mkdir()
        (lock_dir / "owner.json").write_text(
            json.dumps({"token": "parent-token"}) + "\n",
            encoding="utf-8",
        )
        processor = AlignmentProcessor(
            self.input_dir,
            self.output_dir,
            cache_enabled=False,
            acquire_run_lock=False,
            parent_lock_token="parent-token",
        )
        processor.group_sequences_by_multi_id = mock.Mock(return_value={})

        self.assertIsNone(processor.run())
        self.assertTrue(lock_dir.is_dir())

        processor.parent_lock_token = "wrong-token"
        with self.assertRaisesRegex(RuntimeError, "does not match"):
            processor.run()

        (lock_dir / "owner.json").unlink()
        lock_dir.rmdir()

    def test_symlink_aliases_share_one_canonical_workspace_lock(self):
        real_workspace = self.temp_path / "real-workspace"
        real_input = real_workspace / "input"
        real_output = real_workspace / "output"
        real_cache = real_workspace / "cache"
        real_input.mkdir(parents=True)
        alias_workspace = self.temp_path / "workspace-alias"
        alias_workspace.symlink_to(real_workspace, target_is_directory=True)

        canonical_processor = AlignmentProcessor(
            real_input,
            real_output,
            cache_dir=real_cache,
            cache_enabled=False,
        )
        alias_processor = AlignmentProcessor(
            alias_workspace / "input",
            alias_workspace / "output",
            cache_dir=alias_workspace / "cache",
            cache_enabled=False,
        )

        self.assertEqual(
            alias_processor.output_dir,
            canonical_processor.output_dir,
        )
        self.assertEqual(
            alias_processor.run_lock_path,
            canonical_processor.run_lock_path,
        )
        self.assertEqual(
            alias_processor.cache_dir,
            canonical_processor.cache_dir,
        )

    def test_standard_group_file_rejects_foreign_record(self):
        (self.input_dir / "Multi_1.fa").write_text(
            ">Multi_2_transcript_1\nGGCC\n",
            encoding="utf-8",
        )
        groups = self.processor.group_sequences_by_multi_id()

        with self.assertRaisesRegex(RuntimeError, "expected Multi_1"):
            self.processor._load_group_sequences(
                "Multi_1", groups["Multi_1"]
            )

    def test_empty_input_sequence_is_rejected_before_publication(self):
        (self.input_dir / "Multi_1.fa").write_text(
            ">Multi_1_transcript_1\n",
            encoding="utf-8",
        )
        previous_output = (
            self.output_dir / "Multi_previous_aligned_consensus.fa"
        )
        previous_output.write_text(
            ">Multi_previous\nACGT\n", encoding="utf-8"
        )

        with self.assertRaisesRegex(
            RuntimeError,
            "Empty input sequence.*Multi_1",
        ):
            self.processor.run()

        self.assertEqual(
            previous_output.read_text(encoding="utf-8"),
            ">Multi_previous\nACGT\n",
        )

    def test_public_discovery_output_can_be_processed_directly(self):
        (self.input_dir / "Multi_1.fa").write_text(
            ">Multi_1_transcript_1\nACGT\n",
            encoding="utf-8",
        )
        groups = self.processor.group_sequences_by_multi_id()
        self.processor._compute_consensus_record = mock.Mock(
            return_value=self._consensus("Multi_1", "ACGT")
        )

        output_path = self.processor.process_multi_gene_group(
            "Multi_1", groups["Multi_1"]
        )

        self.assertTrue(output_path.is_file())
        self.assertIn(
            ">Multi_1 test consensus",
            output_path.read_text(encoding="utf-8"),
        )

    def test_successful_zero_group_run_removes_stale_outputs(self):
        stale_output = (
            self.output_dir / "Multi_old_aligned_consensus.fa"
        )
        stale_output.write_text(
            ">Multi_old\nAAAA\n", encoding="utf-8"
        )

        result = self.processor.run()

        self.assertIsNone(result)
        self.assertFalse(stale_output.exists())

    def test_interrupted_generation_publication_restores_previous_set(self):
        previous_outputs = {
            "Multi_1": "OLD_ONE",
            "Multi_2": "OLD_TWO",
        }
        for multi_id, sequence in previous_outputs.items():
            (
                self.output_dir / f"{multi_id}_aligned_consensus.fa"
            ).write_text(
                f">{multi_id}\n{sequence}\n", encoding="utf-8"
            )
        staging_dir = self.output_dir / ".test-staging"
        staging_dir.mkdir()
        for multi_id in previous_outputs:
            (
                staging_dir / f"{multi_id}_aligned_consensus.fa"
            ).write_text(
                f">{multi_id}\nNEW\n", encoding="utf-8"
            )

        real_replace = os.replace
        new_publication_count = 0

        def interrupt_second_new_output(source, destination):
            nonlocal new_publication_count
            source = Path(source)
            destination = Path(destination)
            if (
                source.parent == staging_dir
                and destination.parent == self.output_dir
            ):
                new_publication_count += 1
                if new_publication_count == 2:
                    raise KeyboardInterrupt()
            return real_replace(source, destination)

        with mock.patch(
            "src.docker.run_alignment_and_mask.run_alignment_and_mask.os.replace",
            side_effect=interrupt_second_new_output,
        ):
            with self.assertRaises(KeyboardInterrupt):
                self.processor._publish_output_generation(
                    staging_dir, self.output_dir
                )

        for multi_id, sequence in previous_outputs.items():
            output_path = (
                self.output_dir / f"{multi_id}_aligned_consensus.fa"
            )
            self.assertEqual(
                output_path.read_text(encoding="utf-8"),
                f">{multi_id}\n{sequence}\n",
            )
        self.assertEqual(
            list(self.output_dir.glob(
                ".previous-alignment-generation-*"
            )),
            [],
        )

    def test_standalone_processor_uses_shared_workspace_lock(self):
        self.processor.group_sequences_by_multi_id = mock.Mock(
            return_value={}
        )
        completed = threading.Event()
        failures = []

        def execute():
            try:
                self.processor.run()
            except Exception as error:
                failures.append(error)
            finally:
                completed.set()

        with annotation_workspace_lock(self.temp_path):
            worker = threading.Thread(target=execute)
            worker.start()
            time.sleep(0.1)
            self.assertFalse(completed.is_set())

        worker.join(timeout=3)
        self.assertFalse(worker.is_alive())
        self.assertEqual(failures, [])
        self.assertTrue(completed.is_set())

    def test_cache_reuse_restores_independent_atomic_output(self):
        self._set_groups({"Multi_1": self.sequences})
        compute_mock = mock.Mock(
            return_value=self._consensus("Multi_1")
        )
        self.processor._compute_consensus_record = compute_mock

        self.processor.run()
        cache_key = self.processor._group_cache_key(
            "Multi_1", self.sequences
        )
        cache_object, metadata_path, _ = self.processor._cache_paths(
            cache_key
        )
        final_output = (
            self.output_dir / "Multi_1_aligned_consensus.fa"
        )
        self.assertTrue(cache_object.is_file())
        self.assertTrue(metadata_path.is_file())
        self.assertNotEqual(cache_object.stat().st_ino, final_output.stat().st_ino)

        final_output.write_text("damaged final output\n", encoding="utf-8")
        self.processor.run()

        self.assertEqual(compute_mock.call_count, 1)
        self.assertEqual(final_output.read_bytes(), cache_object.read_bytes())
        self.assertNotEqual(cache_object.stat().st_ino, final_output.stat().st_ino)

    def test_cache_key_invalidates_on_order_sequence_and_tool_identity(self):
        baseline_key = self.processor._group_cache_key(
            "Multi_1", self.sequences
        )
        reversed_key = self.processor._group_cache_key(
            "Multi_1", list(reversed(self.sequences))
        )
        changed_sequences = [
            SeqRecord(Seq("ACGA"), id="Multi_1_transcript_1"),
            self.sequences[1],
        ]
        changed_sequence_key = self.processor._group_cache_key(
            "Multi_1", changed_sequences
        )
        other_tool = AlignmentProcessor(
            self.input_dir,
            self.temp_path / "other_output",
            cache_dir=self.cache_dir,
            tool_identity="other-image@sha256:1234",
        )
        changed_tool_key = other_tool._group_cache_key(
            "Multi_1", self.sequences
        )

        self.assertNotEqual(baseline_key, reversed_key)
        self.assertNotEqual(baseline_key, changed_sequence_key)
        self.assertNotEqual(baseline_key, changed_tool_key)

    def test_corrupt_or_missing_cache_files_are_recomputed(self):
        self._set_groups({"Multi_1": self.sequences})
        compute_mock = mock.Mock(
            return_value=self._consensus("Multi_1")
        )
        self.processor._compute_consensus_record = compute_mock

        self.processor.run()
        cache_key = self.processor._group_cache_key(
            "Multi_1", self.sequences
        )
        cache_object, metadata_path, _ = self.processor._cache_paths(
            cache_key
        )

        cache_object.write_text(
            ">Multi_1\nCORRUPT\n", encoding="utf-8"
        )
        self.processor.run()
        self.assertEqual(compute_mock.call_count, 2)

        metadata_path.unlink()
        self.processor.run()
        self.assertEqual(compute_mock.call_count, 3)
        with open(metadata_path, "r", encoding="utf-8") as handle:
            metadata = json.load(handle)
        self.assertEqual(
            metadata["output_sha256"],
            self.processor._file_sha256(cache_object),
        )

    def test_bounded_parallel_run_processes_groups_concurrently(self):
        processor = AlignmentProcessor(
            self.input_dir,
            self.output_dir,
            cache_dir=self.cache_dir,
            workers=2,
            cache_enabled=False,
        )
        groups = {
            f"Multi_{index}": [
                SeqRecord(
                    Seq("ACGT"), id=f"Multi_{index}_transcript_1"
                )
            ]
            for index in range(1, 5)
        }
        processor.group_sequences_by_multi_id = mock.Mock(return_value=groups)

        state_lock = threading.Lock()
        active_count = 0
        maximum_active = 0
        first_pair = threading.Barrier(2)

        def process_group(multi_id, sequences):
            nonlocal active_count, maximum_active
            with state_lock:
                active_count += 1
                maximum_active = max(maximum_active, active_count)
            try:
                if multi_id in {"Multi_1", "Multi_2"}:
                    first_pair.wait(timeout=2)
                time.sleep(0.02)
                output_path = (
                    processor.output_dir
                    / f"{multi_id}_aligned_consensus.fa"
                )
                output_path.write_text(
                    f">{multi_id}\nACGT\n", encoding="utf-8"
                )
                return output_path
            finally:
                with state_lock:
                    active_count -= 1

        processor.process_multi_gene_group = mock.Mock(
            side_effect=process_group
        )

        result = processor.run()

        self.assertEqual(result["successful_groups"], 4)
        self.assertEqual(processor.process_multi_gene_group.call_count, 4)
        self.assertEqual(maximum_active, 2)

    def test_run_aggregates_failures_and_keeps_successful_cache(self):
        groups = {
            "Multi_1": self.sequences,
            "Multi_2": [
                SeqRecord(Seq("GGCC"), id="Multi_2_transcript_1")
            ],
        }
        self._set_groups(groups)
        attempted_groups = []

        def compute_group(multi_id, sequences):
            attempted_groups.append(multi_id)
            if multi_id == "Multi_2":
                raise RuntimeError("simulated MAFFT failure")
            return self._consensus(multi_id)

        compute_mock = mock.Mock(side_effect=compute_group)
        self.processor._compute_consensus_record = compute_mock
        previous_output = (
            self.output_dir / "Multi_previous_aligned_consensus.fa"
        )
        previous_output.write_text(
            ">Multi_previous\nTTTT\n", encoding="utf-8"
        )

        with self.assertRaisesRegex(RuntimeError, "1 of 2 groups"):
            self.processor.run()

        self.assertEqual(set(attempted_groups), {"Multi_1", "Multi_2"})
        self.assertEqual(
            previous_output.read_text(encoding="utf-8"),
            ">Multi_previous\nTTTT\n",
        )
        self.processor.group_sequences_by_multi_id.return_value = {
            "Multi_1": self.sequences
        }
        compute_mock.reset_mock()
        self.processor.run()
        compute_mock.assert_not_called()

    def test_invalid_workers_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "positive integer"):
            AlignmentProcessor(
                self.input_dir,
                self.output_dir,
                workers=0,
            )


if __name__ == "__main__":
    unittest.main()
