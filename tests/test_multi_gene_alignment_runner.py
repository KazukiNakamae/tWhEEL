import os
import subprocess
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock

from src.modules.multi_gene_alignment import (
    MultiGeneAlignmentRunner,
    WorkspaceLockSafetyError,
)


class CommandHarness:
    """Simulate Docker image and container commands for runner tests."""

    def __init__(
        self,
        output_dir,
        output_ids=(),
        image_exists=True,
        image_identity='sha256:existing-image',
        empty_identity=False,
    ):
        self.output_dir = Path(output_dir)
        self.output_ids = tuple(output_ids)
        self.image_exists = image_exists
        self.image_identity = image_identity
        self.empty_identity = empty_identity
        self.calls = []

    def __call__(self, command, **kwargs):
        command = list(command)
        self.calls.append((command, kwargs))

        if command[:3] == ['docker', 'image', 'inspect']:
            if not self.image_exists:
                return subprocess.CompletedProcess(
                    command,
                    1,
                    stdout='',
                    stderr='No such image',
                )
            stdout = '' if self.empty_identity else self.image_identity + '\n'
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=stdout,
                stderr='',
            )

        if command[:2] == ['docker', 'build']:
            self.image_exists = True
            return subprocess.CompletedProcess(
                command,
                0,
                stdout='Successfully built\n',
                stderr='',
            )

        if command[:2] == ['docker', 'run']:
            self.output_dir.mkdir(parents=True, exist_ok=True)
            for multi_id in self.output_ids:
                output_path = (
                    self.output_dir
                    / f'{multi_id}_aligned_consensus.fa'
                )
                output_path.write_text(
                    f'>{multi_id}\nACGT\n',
                    encoding='utf-8',
                )
            return subprocess.CompletedProcess(
                command,
                0,
                stdout='processed=1 cached=0 failed=0\n',
                stderr='',
            )

        if command[:3] == ['docker', 'rm', '-f']:
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=command[-1] + '\n',
                stderr='',
            )

        raise AssertionError(f'Unexpected command: {command}')

    def commands_starting_with(self, prefix):
        return [
            command
            for command, _ in self.calls
            if command[:len(prefix)] == list(prefix)
        ]


class MultiGeneAlignmentRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_dir.cleanup)
        self.root = Path(self.temp_dir.name)
        self.intermediate_dir = self.root / 'intermediate'
        self.input_dir = self.intermediate_dir / 'multi_genes'
        self.output_dir = self.intermediate_dir / 'multi_genes_aligned'
        self.input_dir.mkdir(parents=True)
        self.build_context = self.root / 'build-context'
        self.build_context.mkdir()
        (self.build_context / 'Dockerfile').write_text(
            'FROM scratch\n', encoding='utf-8'
        )
        (self.build_context / 'requirements.txt').write_text(
            'biopython\n', encoding='utf-8'
        )
        (self.build_context / 'runner.py').write_text(
            'print("runner")\n', encoding='utf-8'
        )
        self.logger = mock.Mock()

    def make_runner(self, harness, settings=None):
        config = {}
        if settings is not None:
            config['multi_gene_alignment'] = settings
        runner = MultiGeneAlignmentRunner(
            self.intermediate_dir,
            config,
            self.logger,
            run_command=harness,
        )
        runner.build_context = self.build_context
        return runner

    def test_existing_image_skips_build_and_uses_actual_image_identity(self):
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
        )
        runner = self.make_runner(harness)

        outputs = runner.run({'Multi_1'})

        self.assertEqual(
            outputs,
            [
                str(
                    (
                        self.output_dir
                        / 'Multi_1_aligned_consensus.fa'
                    ).resolve()
                )
            ],
        )
        self.assertEqual(
            harness.commands_starting_with(['docker', 'build']),
            [],
        )
        run_command = harness.commands_starting_with(['docker', 'run'])[0]
        self.assertEqual(
            run_command[run_command.index('--workers') + 1],
            '4',
        )
        self.assertEqual(
            run_command[run_command.index('--tool-identity') + 1],
            'sha256:existing-image',
        )
        image_reference_index = run_command.index('/data') + 1
        self.assertEqual(
            run_command[image_reference_index],
            'sha256:existing-image',
        )
        self.assertEqual(run_command[run_command.index('-m') + 1], '4g')
        self.assertEqual(
            run_command[run_command.index('--cache-dir') + 1],
            'multi_genes_alignment_cache',
        )
        self.assertIn('--no-run-lock', run_command)
        parent_lock_token = run_command[
            run_command.index('--parent-lock-token') + 1
        ]
        self.assertRegex(parent_lock_token, r'^[0-9a-f]{32}$')
        self.assertTrue(runner.cache_dir.is_dir())
        self.assertEqual(
            list(self.intermediate_dir.glob(
                '.multi-gene-alignment-backup-*'
            )),
            [],
        )
        self.assertEqual(
            len(harness.commands_starting_with(['docker', 'rm', '-f'])),
            0,
        )
        self.assertFalse(
            (
                self.intermediate_dir / '.annotation-workspace.lock.d'
            ).exists()
        )
        self.assertTrue(any(
            call.args
            and call.args[0] == 'Multi-gene alignment summary:\n%s'
            for call in self.logger.info.call_args_list
        ))

    def test_unconfirmed_container_cleanup_preserves_workspace_lock(self):
        self.output_dir.mkdir(parents=True)
        previous_output = self.output_dir / 'Multi_1_aligned_consensus.fa'
        previous_output.write_text(
            '>Multi_1\nPREVIOUS\n', encoding='utf-8'
        )
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
        )

        def unsafe_cleanup_harness(command, **kwargs):
            command = list(command)
            if command[:2] == ['docker', 'run']:
                raise KeyboardInterrupt()
            if command[:3] == ['docker', 'rm', '-f']:
                return subprocess.CompletedProcess(
                    command,
                    1,
                    stdout='',
                    stderr='Cannot connect to the Docker daemon',
                )
            return harness(command, **kwargs)

        runner = MultiGeneAlignmentRunner(
            self.intermediate_dir,
            {},
            self.logger,
            run_command=unsafe_cleanup_harness,
        )
        runner.build_context = self.build_context

        with self.assertRaises(WorkspaceLockSafetyError):
            runner.run({'Multi_1'})

        self.assertFalse(previous_output.exists())
        backup_directories = list(self.intermediate_dir.glob(
            '.multi-gene-alignment-backup-*'
        ))
        self.assertEqual(len(backup_directories), 1)
        backup_output = (
            backup_directories[0]
            / 'Multi_1_aligned_consensus.fa'
        )
        self.assertEqual(
            backup_output.read_text(encoding='utf-8'),
            '>Multi_1\nPREVIOUS\n',
        )
        previous_output.write_text(
            '>Multi_1\nORPHAN_NEW\n', encoding='utf-8'
        )
        self.assertEqual(
            backup_output.read_text(encoding='utf-8'),
            '>Multi_1\nPREVIOUS\n',
        )
        lock_dir = (
            self.intermediate_dir / '.annotation-workspace.lock.d'
        )
        self.assertTrue(lock_dir.is_dir())
        self.assertTrue((lock_dir / 'owner.json').is_file())

    def test_missing_image_is_built_once_without_sudo(self):
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
            image_exists=False,
        )
        runner = self.make_runner(harness)

        runner.run({'Multi_1'})

        build_commands = harness.commands_starting_with(['docker', 'build'])
        self.assertEqual(len(build_commands), 1)
        build_command = build_commands[0]
        self.assertNotIn('sudo', build_command)
        self.assertEqual(build_command[-1], str(self.build_context))
        build_call = next(
            kwargs
            for command, kwargs in harness.calls
            if command[:2] == ['docker', 'build']
        )
        self.assertEqual(
            build_call['env']['BUILDX_CONFIG'],
            str(runner.buildx_config_dir),
        )
        self.assertTrue(runner.buildx_config_dir.is_dir())
        image_tag = build_command[build_command.index('-t') + 1]
        self.assertRegex(
            image_tag,
            r'^twheel-run-alignment-and-mask:[0-9a-f]{64}$',
        )
        self.assertEqual(
            len(harness.commands_starting_with([
                'docker', 'image', 'inspect'
            ])),
            3,
        )
        run_command = harness.commands_starting_with(['docker', 'run'])[0]
        self.assertEqual(
            run_command[run_command.index('--tool-identity') + 1],
            'sha256:existing-image',
        )
        image_reference_index = run_command.index('/data') + 1
        self.assertEqual(
            run_command[image_reference_index],
            'sha256:existing-image',
        )

    def test_successful_inspect_requires_a_nonempty_image_identity(self):
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
            empty_identity=True,
        )
        runner = self.make_runner(harness)

        with self.assertRaisesRegex(RuntimeError, 'empty identity'):
            runner.run({'Multi_1'})

        self.assertEqual(
            harness.commands_starting_with(['docker', 'run']),
            [],
        )

    def test_reinspect_under_lock_avoids_duplicate_build(self):
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
            image_exists=False,
        )
        original_call = harness.__call__
        inspect_count = 0

        def image_appears_during_lock(command, **kwargs):
            nonlocal inspect_count
            if list(command)[:3] == ['docker', 'image', 'inspect']:
                inspect_count += 1
                if inspect_count == 2:
                    harness.image_exists = True
            return original_call(command, **kwargs)

        runner = MultiGeneAlignmentRunner(
            self.intermediate_dir,
            {},
            self.logger,
            run_command=image_appears_during_lock,
        )
        runner.build_context = self.build_context

        runner.run({'Multi_1'})

        self.assertEqual(inspect_count, 2)
        self.assertEqual(
            harness.commands_starting_with(['docker', 'build']),
            [],
        )

    def test_custom_settings_and_disabled_cache_reach_docker(self):
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
        )
        runner = self.make_runner(harness, {
            'workers': 7,
            'cache_enabled': False,
            'memory': '512m',
        })

        runner.run({'Multi_1'})

        run_command = harness.commands_starting_with(['docker', 'run'])[0]
        self.assertEqual(
            run_command[run_command.index('--workers') + 1],
            '7',
        )
        self.assertEqual(run_command[run_command.index('-m') + 1], '512m')
        self.assertIn('--no-cache', run_command)
        self.assertNotIn('--cache-dir', run_command)

    def test_none_expected_ids_are_derived_from_input_files(self):
        (self.input_dir / 'Multi_1.fa').write_text(
            '>Multi_1_transcript\nACGT\n',
            encoding='utf-8',
        )
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
        )
        runner = self.make_runner(harness)

        outputs = runner.run(None)

        self.assertEqual(
            outputs,
            [
                str(
                    (
                        self.output_dir
                        / 'Multi_1_aligned_consensus.fa'
                    ).resolve()
                )
            ],
        )

    def test_failed_generation_restores_prior_complete_output_set(self):
        self.output_dir.mkdir(parents=True)
        stale_path = self.output_dir / 'Multi_old_aligned_consensus.fa'
        stale_path.write_text('>stale\nAAAA\n', encoding='utf-8')
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_unexpected'],
        )
        runner = self.make_runner(harness)

        with self.assertRaisesRegex(
            RuntimeError,
            r'missing=Multi_expected; unexpected=Multi_unexpected',
        ):
            runner.run({'Multi_expected'})

        self.assertEqual(
            stale_path.read_text(encoding='utf-8'),
            '>stale\nAAAA\n',
        )
        self.assertFalse(
            (
                self.output_dir
                / 'Multi_unexpected_aligned_consensus.fa'
            ).exists()
        )

    def test_invalid_consensus_content_rolls_back_previous_valid_set(self):
        invalid_contents = (
            ('empty file', ''),
            ('header only', '>Multi_1\n'),
            ('wrong identifier', '>Wrong_ID\nACGT\n'),
            (
                'multiple records',
                '>Multi_1\nACGT\n>Multi_1_extra\nACGT\n',
            ),
        )
        for label, invalid_content in invalid_contents:
            with self.subTest(label=label):
                self.output_dir.mkdir(parents=True, exist_ok=True)
                previous_output = (
                    self.output_dir / 'Multi_1_aligned_consensus.fa'
                )
                previous_output.write_text(
                    '>Multi_1\nPREVIOUS\n', encoding='utf-8'
                )
                harness = CommandHarness(
                    self.output_dir,
                    output_ids=['Multi_1'],
                )

                def invalid_output_harness(command, **kwargs):
                    result = harness(command, **kwargs)
                    if list(command)[:2] == ['docker', 'run']:
                        previous_output.write_text(
                            invalid_content, encoding='utf-8'
                        )
                    return result

                runner = MultiGeneAlignmentRunner(
                    self.intermediate_dir,
                    {},
                    self.logger,
                    run_command=invalid_output_harness,
                )
                runner.build_context = self.build_context

                with self.assertRaisesRegex(
                    RuntimeError,
                    'Invalid Multi-gene consensus FASTA',
                ):
                    runner.run({'Multi_1'})

                self.assertEqual(
                    previous_output.read_text(encoding='utf-8'),
                    '>Multi_1\nPREVIOUS\n',
                )

    def test_keyboard_interrupt_restores_prior_complete_output_set(self):
        self.output_dir.mkdir(parents=True)
        previous_output = self.output_dir / 'Multi_1_aligned_consensus.fa'
        previous_output.write_text(
            '>Multi_1\nPREVIOUS\n', encoding='utf-8'
        )
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
        )

        def interrupting_harness(command, **kwargs):
            if list(command)[:2] == ['docker', 'run']:
                raise KeyboardInterrupt()
            return harness(command, **kwargs)

        runner = MultiGeneAlignmentRunner(
            self.intermediate_dir,
            {},
            self.logger,
            run_command=interrupting_harness,
        )
        runner.build_context = self.build_context

        with self.assertRaises(KeyboardInterrupt):
            runner.run({'Multi_1'})

        self.assertEqual(
            previous_output.read_text(encoding='utf-8'),
            '>Multi_1\nPREVIOUS\n',
        )
        self.assertEqual(
            list(self.intermediate_dir.glob(
                '.multi-gene-alignment-backup-*'
            )),
            [],
        )
        self.assertEqual(
            len(harness.commands_starting_with(['docker', 'rm', '-f'])),
            1,
        )
        self.assertFalse(
            (
                self.intermediate_dir / '.annotation-workspace.lock.d'
            ).exists()
        )

    def test_partial_backup_failure_restores_every_moved_output(self):
        self.output_dir.mkdir(parents=True)
        first_output = self.output_dir / 'Multi_1_aligned_consensus.fa'
        second_output = self.output_dir / 'Multi_2_aligned_consensus.fa'
        first_output.write_text('first\n', encoding='utf-8')
        second_output.write_text('second\n', encoding='utf-8')
        runner = self.make_runner(CommandHarness(self.output_dir))
        real_replace = os.replace
        backup_move_count = 0

        def fail_second_backup_move(source, destination):
            nonlocal backup_move_count
            source = Path(source)
            destination = Path(destination)
            if (
                source.parent.resolve() == self.output_dir.resolve()
                and destination.parent.name.startswith(
                    '.multi-gene-alignment-backup-'
                )
            ):
                backup_move_count += 1
                if backup_move_count == 2:
                    raise OSError('simulated backup move failure')
            return real_replace(source, destination)

        with mock.patch(
            'src.modules.multi_gene_alignment.os.replace',
            side_effect=fail_second_backup_move,
        ):
            with self.assertRaisesRegex(OSError, 'backup move failure'):
                runner._backup_published_outputs()

        self.assertEqual(first_output.read_text(encoding='utf-8'), 'first\n')
        self.assertEqual(second_output.read_text(encoding='utf-8'), 'second\n')
        self.assertEqual(
            list(self.intermediate_dir.glob(
                '.multi-gene-alignment-backup-*'
            )),
            [],
        )

    def test_failed_restore_preserves_backup_for_manual_recovery(self):
        self.output_dir.mkdir(parents=True)
        previous_output = self.output_dir / 'Multi_1_aligned_consensus.fa'
        previous_output.write_text(
            '>Multi_1\nPREVIOUS\n', encoding='utf-8'
        )
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
        )

        def failing_harness(command, **kwargs):
            if list(command)[:2] == ['docker', 'run']:
                raise RuntimeError('simulated runner failure')
            return harness(command, **kwargs)

        runner = MultiGeneAlignmentRunner(
            self.intermediate_dir,
            {},
            self.logger,
            run_command=failing_harness,
        )
        runner.build_context = self.build_context

        with mock.patch.object(
            runner,
            '_restore_published_outputs',
            side_effect=OSError('simulated restore failure'),
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                'Preserved backup:',
            ):
                runner.run({'Multi_1'})

        backup_directories = list(self.intermediate_dir.glob(
            '.multi-gene-alignment-backup-*'
        ))
        self.assertEqual(len(backup_directories), 1)
        self.assertEqual(
            (
                backup_directories[0]
                / 'Multi_1_aligned_consensus.fa'
            ).read_text(encoding='utf-8'),
            '>Multi_1\nPREVIOUS\n',
        )
        self.assertTrue(
            (
                self.intermediate_dir / '.annotation-workspace.lock.d'
            ).is_dir()
        )

    def test_same_intermediate_directory_serializes_complete_docker_runs(self):
        harness = CommandHarness(
            self.output_dir,
            output_ids=['Multi_1'],
        )
        first_run_entered = threading.Event()
        release_first_run = threading.Event()
        call_lock = threading.Lock()
        docker_run_count = 0

        def blocking_harness(command, **kwargs):
            nonlocal docker_run_count
            if list(command)[:2] == ['docker', 'run']:
                with call_lock:
                    docker_run_count += 1
                    current_run = docker_run_count
                if current_run == 1:
                    first_run_entered.set()
                    if not release_first_run.wait(timeout=3):
                        raise RuntimeError('test timed out releasing first run')
            return harness(command, **kwargs)

        first_runner = MultiGeneAlignmentRunner(
            self.intermediate_dir,
            {},
            self.logger,
            run_command=blocking_harness,
        )
        second_runner = MultiGeneAlignmentRunner(
            self.intermediate_dir,
            {},
            self.logger,
            run_command=blocking_harness,
        )
        first_runner.build_context = self.build_context
        second_runner.build_context = self.build_context
        failures = []

        def run_alignment(runner):
            try:
                runner.run({'Multi_1'})
            except Exception as error:
                failures.append(error)

        first_thread = threading.Thread(
            target=run_alignment, args=(first_runner,)
        )
        second_thread = threading.Thread(
            target=run_alignment, args=(second_runner,)
        )
        first_thread.start()
        self.assertTrue(first_run_entered.wait(timeout=3))
        second_thread.start()
        time.sleep(0.1)
        self.assertEqual(docker_run_count, 1)
        release_first_run.set()
        first_thread.join(timeout=3)
        second_thread.join(timeout=3)

        self.assertFalse(first_thread.is_alive())
        self.assertFalse(second_thread.is_alive())
        self.assertEqual(failures, [])
        self.assertEqual(docker_run_count, 2)

    def test_build_context_hash_uses_sources_and_ignores_python_cache(self):
        harness = CommandHarness(self.output_dir)
        runner = self.make_runner(harness)
        first_digest = runner._build_context_digest()

        pycache = self.build_context / '__pycache__'
        pycache.mkdir()
        (pycache / 'runner.cpython-312.pyc').write_bytes(b'ignored bytecode')
        (self.build_context / '.DS_Store').write_bytes(b'ignored metadata')
        self.assertEqual(runner._build_context_digest(), first_digest)

        (self.build_context / 'runner.py').write_text(
            'print("updated")\n', encoding='utf-8'
        )
        self.assertNotEqual(runner._build_context_digest(), first_digest)

    def test_invalid_configuration_is_rejected(self):
        invalid_sections = [
            [],
            {'workers': 0},
            {'workers': True},
            {'workers': '4'},
            {'cache_enabled': 1},
            {'cache_enabled': 'true'},
            {'memory': ''},
            {'memory': '0g'},
            {'memory': '4'},
            {'memory': 4},
        ]

        for section in invalid_sections:
            with self.subTest(section=section):
                with self.assertRaises((TypeError, ValueError)):
                    MultiGeneAlignmentRunner(
                        self.intermediate_dir,
                        {'multi_gene_alignment': section},
                        self.logger,
                        run_command=lambda *args, **kwargs: None,
                    )

        with self.assertRaises(TypeError):
            MultiGeneAlignmentRunner(
                self.intermediate_dir,
                None,
                self.logger,
                run_command=lambda *args, **kwargs: None,
            )


if __name__ == '__main__':
    unittest.main()
