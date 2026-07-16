import multiprocessing
import subprocess
import tempfile
import threading
import time
import unittest
from pathlib import Path
from unittest import mock

from src.modules.annotation import Annotation
from src.modules.multi_gene_alignment import annotation_workspace_lock
from src.modules.multi_gene_alignment import WorkspaceLockSafetyError
from src.modules.sgd_annotation import Sgd_annotation


def _hold_workspace_lock(intermediate_dir, entered, release):
    """Test helper that retains the shared lock in a child process."""
    with annotation_workspace_lock(intermediate_dir):
        entered.set()
        if not release.wait(timeout=5):
            raise RuntimeError('workspace-lock release timed out')


def _enter_workspace_lock(intermediate_dir, entered):
    """Test helper that reports entry into the shared lock."""
    with annotation_workspace_lock(intermediate_dir):
        entered.set()


class MultiGeneAlignmentContractTests(unittest.TestCase):
    def test_workspace_generation_lock_is_reentrant_in_one_thread(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            entered_nested_scope = False
            with annotation_workspace_lock(temp_dir):
                with annotation_workspace_lock(temp_dir):
                    entered_nested_scope = True
            self.assertTrue(entered_nested_scope)

    def test_nested_safety_error_preserves_outer_workspace_lock(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            with annotation_workspace_lock(temp_dir):
                try:
                    with annotation_workspace_lock(temp_dir):
                        raise WorkspaceLockSafetyError('unsafe container state')
                except WorkspaceLockSafetyError:
                    pass
            lock_dir = (
                Path(temp_dir) / '.annotation-workspace.lock.d'
            )
            self.assertTrue(lock_dir.is_dir())
            self.assertTrue((lock_dir / 'owner.json').is_file())

    def test_workspace_generation_lock_serializes_os_processes(self):
        context = multiprocessing.get_context('fork')
        with tempfile.TemporaryDirectory() as temp_dir:
            first_entered = context.Event()
            second_entered = context.Event()
            release_first = context.Event()
            first_process = context.Process(
                target=_hold_workspace_lock,
                args=(temp_dir, first_entered, release_first),
            )
            second_process = context.Process(
                target=_enter_workspace_lock,
                args=(temp_dir, second_entered),
            )
            first_process.start()
            self.assertTrue(first_entered.wait(timeout=3))
            second_process.start()
            try:
                self.assertFalse(second_entered.wait(timeout=0.2))
                release_first.set()
                self.assertTrue(second_entered.wait(timeout=3))
            finally:
                release_first.set()
                first_process.join(timeout=3)
                second_process.join(timeout=3)
                if first_process.is_alive():
                    first_process.terminate()
                    first_process.join(timeout=3)
                if second_process.is_alive():
                    second_process.terminate()
                    second_process.join(timeout=3)
            self.assertEqual(first_process.exitcode, 0)
            self.assertEqual(second_process.exitcode, 0)

    def test_forked_child_does_not_inherit_reentrant_lock_depth(self):
        context = multiprocessing.get_context('fork')
        with tempfile.TemporaryDirectory() as temp_dir:
            child_entered = context.Event()
            child_process = context.Process(
                target=_enter_workspace_lock,
                args=(temp_dir, child_entered),
            )
            with annotation_workspace_lock(temp_dir):
                child_process.start()
                self.assertFalse(child_entered.wait(timeout=0.2))
            try:
                self.assertTrue(child_entered.wait(timeout=3))
            finally:
                child_process.join(timeout=3)
                if child_process.is_alive():
                    child_process.terminate()
                    child_process.join(timeout=3)
            self.assertEqual(child_process.exitcode, 0)

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
            def successful_without_alignment_output(command, **kwargs):
                command = list(command)
                stdout = (
                    'sha256:test-image\n'
                    if command[:3] == ['docker', 'image', 'inspect']
                    else ''
                )
                return subprocess.CompletedProcess(
                    args=command,
                    returncode=0,
                    stdout=stdout,
                    stderr='',
                )

            with mock.patch(
                run_target,
                side_effect=successful_without_alignment_output,
            ):
                with self.assertRaisesRegex(
                    RuntimeError, 'missing=Multi_1'
                ):
                    module._run_alignment_process({'Multi_1'})
            self.assertEqual(
                stale_output.read_text(encoding='utf-8'),
                '>Multi_1\nSTALE\n',
            )

    def test_empty_categories_remove_stale_publications_but_keep_cache(self):
        header = (
            'query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\t'
            'description\tgene_name\tclassification\tunique_id\n'
        )
        single_only = (
            header
            + 'tx1.p1\ts1\t99\t4\t0\t10\tdesc\tGENE\tSingle\t'
            'Single_1\n'
        )
        multi_only = (
            header
            + 'tx1.p1\ts1\t99\t4\t0\t10\tdesc\tGENE\tMulti\t'
            'Multi_1\n'
            + 'tx2.p1\ts1\t99\t4\t0\t10\tdesc\tGENE\tMulti\t'
            'Multi_1\n'
        )

        for module_class in (Annotation, Sgd_annotation):
            with self.subTest(module_class=module_class.__name__):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    intermediate = Path(temp_dir_name)
                    (intermediate / 'transcriptome.fa').write_text(
                        '>tx1\nACGT\n>tx2\nGGCC\n',
                        encoding='utf-8',
                    )
                    classified = intermediate / 'classified.tsv'
                    stale_directories = (
                        intermediate / 'multi_genes',
                        intermediate / 'multi_genes_aligned',
                        intermediate / 'multi_genes_aligned_oneline',
                    )
                    for directory in stale_directories:
                        directory.mkdir()
                    (intermediate / 'multi_genes.fa').write_text(
                        '>Multi_old\nAAAA\n', encoding='utf-8'
                    )
                    (intermediate / 'multi_genes' / 'Multi_old.fa').write_text(
                        '>Multi_old_tx\nAAAA\n', encoding='utf-8'
                    )
                    (
                        intermediate
                        / 'multi_genes_aligned'
                        / 'Multi_old_aligned_consensus.fa'
                    ).write_text('>Multi_old\nAAAA\n', encoding='utf-8')
                    (
                        intermediate
                        / 'multi_genes_aligned_oneline'
                        / 'Multi_old_aligned_consensus.fa'
                    ).write_text('>Multi_old\nAAAA\n', encoding='utf-8')
                    (intermediate / 'multi_genes_aligned_oneline.fa').write_text(
                        '>Multi_old\nAAAA\n', encoding='utf-8'
                    )
                    cache_file = (
                        intermediate
                        / 'multi_genes_alignment_cache'
                        / 'objects'
                        / 'cache.fa'
                    )
                    cache_file.parent.mkdir(parents=True)
                    cache_file.write_text('cached\n', encoding='utf-8')

                    classified.write_text(single_only, encoding='utf-8')
                    module = module_class({
                        'output': {'intermediate': temp_dir_name},
                    })
                    module.extract_gene_sequences(classified)

                    self.assertFalse(
                        (intermediate / 'multi_genes.fa').exists()
                    )
                    for directory in stale_directories:
                        self.assertEqual(list(directory.glob('*.fa')), [])
                    self.assertFalse(
                        (
                            intermediate / 'multi_genes_aligned_oneline.fa'
                        ).exists()
                    )
                    self.assertTrue(cache_file.exists())
                    self.assertTrue(
                        (intermediate / 'single_genes.fa').exists()
                    )

                    classified.write_text(multi_only, encoding='utf-8')
                    with mock.patch.object(
                        module,
                        '_run_alignment_process',
                        return_value=[],
                    ):
                        module.extract_gene_sequences(classified)
                    self.assertFalse(
                        (intermediate / 'single_genes.fa').exists()
                    )

    def test_annotation_run_serializes_shared_workspace_generation(self):
        for module_class in (Annotation, Sgd_annotation):
            with self.subTest(module_class=module_class.__name__):
                with tempfile.TemporaryDirectory() as temp_dir_name:
                    config = {
                        'output': {'intermediate': temp_dir_name},
                    }
                    first_module = module_class(config)
                    second_module = module_class(config)
                    first_entered = threading.Event()
                    second_entered = threading.Event()
                    release_first = threading.Event()
                    failures = []

                    def first_generation():
                        first_entered.set()
                        if not release_first.wait(timeout=3):
                            raise RuntimeError('test release timed out')
                        return 'first'

                    def second_generation():
                        second_entered.set()
                        return 'second'

                    first_module._run_locked = first_generation
                    second_module._run_locked = second_generation

                    def execute(module):
                        try:
                            module.run()
                        except Exception as error:
                            failures.append(error)

                    first_thread = threading.Thread(
                        target=execute, args=(first_module,)
                    )
                    second_thread = threading.Thread(
                        target=execute, args=(second_module,)
                    )
                    first_thread.start()
                    self.assertTrue(first_entered.wait(timeout=3))
                    second_thread.start()
                    time.sleep(0.1)
                    self.assertFalse(second_entered.is_set())
                    release_first.set()
                    first_thread.join(timeout=3)
                    second_thread.join(timeout=3)

                    self.assertFalse(first_thread.is_alive())
                    self.assertFalse(second_thread.is_alive())
                    self.assertEqual(failures, [])
                    self.assertTrue(second_entered.is_set())


if __name__ == '__main__':
    unittest.main()
