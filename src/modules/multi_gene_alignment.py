"""Build and run the cached multi-gene alignment container safely."""

import fcntl
import hashlib
import json
import logging
import os
import re
import shutil
import socket
import subprocess
import tempfile
import threading
import time
import uuid
from collections.abc import Mapping
from contextlib import contextmanager
from pathlib import Path


_WORKSPACE_LOCKS = {}
_WORKSPACE_LOCKS_GUARD = threading.Lock()
_WORKSPACE_LOCK_DEPTHS = threading.local()
_WORKSPACE_LOCK_TIMEOUT_SECONDS = 30.0


class WorkspaceLockSafetyError(RuntimeError):
    """Signal that workspace ownership must remain fail-closed."""

    preserve_workspace_lock = True


def _reset_workspace_lock_state_after_fork():
    """Discard inherited process-local lock bookkeeping in a forked child."""
    global _WORKSPACE_LOCKS
    global _WORKSPACE_LOCKS_GUARD
    global _WORKSPACE_LOCK_DEPTHS
    _WORKSPACE_LOCKS = {}
    _WORKSPACE_LOCKS_GUARD = threading.Lock()
    _WORKSPACE_LOCK_DEPTHS = threading.local()


if hasattr(os, 'register_at_fork'):
    os.register_at_fork(after_in_child=_reset_workspace_lock_state_after_fork)


def _workspace_thread_lock(lock_path):
    """Return a process-local reentrant lock for one workspace path."""
    lock_key = str(Path(lock_path).resolve())
    with _WORKSPACE_LOCKS_GUARD:
        return _WORKSPACE_LOCKS.setdefault(lock_key, threading.RLock())


def _acquire_workspace_directory_lock(lock_dir):
    """Acquire a cross-kernel lock using atomic directory creation."""
    lock_dir = Path(lock_dir)
    lock_dir.parent.mkdir(parents=True, exist_ok=True)
    deadline = time.monotonic() + _WORKSPACE_LOCK_TIMEOUT_SECONDS
    while True:
        try:
            lock_dir.mkdir()
            break
        except FileExistsError:
            if time.monotonic() >= deadline:
                owner_path = lock_dir / 'owner.json'
                try:
                    owner_details = owner_path.read_text(encoding='utf-8').strip()
                except OSError:
                    owner_details = 'owner metadata unavailable'
                raise RuntimeError(
                    'Timed out waiting for the tWhEEL workspace lock. '
                    f'Lock directory: {lock_dir}. Owner: {owner_details}. '
                    'If the owning process or container crashed, verify that '
                    'no tWhEEL container is still running before removing the '
                    'lock directory manually.'
                )
            time.sleep(0.1)

    ownership_token = uuid.uuid4().hex
    owner_path = lock_dir / 'owner.json'
    temporary_owner_path = lock_dir / f'.owner-{ownership_token}.tmp'
    owner_payload = {
        'token': ownership_token,
        'pid': os.getpid(),
        'hostname': socket.gethostname(),
        'created_at_epoch': time.time(),
    }
    try:
        with temporary_owner_path.open('w', encoding='utf-8') as handle:
            json.dump(owner_payload, handle, sort_keys=True)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_owner_path, owner_path)
    except BaseException:
        try:
            temporary_owner_path.unlink()
        except FileNotFoundError:
            pass
        shutil.rmtree(lock_dir, ignore_errors=True)
        raise
    return ownership_token


def _release_workspace_directory_lock(lock_dir, ownership_token):
    """Release only the directory lock owned by this context."""
    lock_dir = Path(lock_dir)
    owner_path = lock_dir / 'owner.json'
    try:
        owner_payload = json.loads(owner_path.read_text(encoding='utf-8'))
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as error:
        raise RuntimeError(
            f'Cannot verify ownership of workspace lock: {lock_dir}'
        ) from error
    if owner_payload.get('token') != ownership_token:
        raise RuntimeError(
            f'Workspace lock ownership changed unexpectedly: {lock_dir}'
        )
    owner_path.unlink()
    lock_dir.rmdir()


@contextmanager
def _exclusive_file_lock(lock_path):
    """Hold an exclusive cross-process lock for the requested scope."""
    lock_path = Path(lock_path)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open('a+b') as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


@contextmanager
def annotation_workspace_lock(intermediate_dir):
    """Serialize annotation generations sharing one intermediate directory."""
    lock_path = (
        Path(intermediate_dir).expanduser().resolve()
        / '.annotation-workspace.lock.d'
    )
    lock_key = str(lock_path)
    thread_lock = _workspace_thread_lock(lock_path)
    with thread_lock:
        depths = getattr(_WORKSPACE_LOCK_DEPTHS, 'values', None)
        if depths is None:
            depths = {}
            _WORKSPACE_LOCK_DEPTHS.values = depths

        lock_state = depths.get(lock_key)
        if lock_state is not None:
            lock_state['depth'] += 1
            try:
                yield lock_state['token']
            except BaseException as error:
                if getattr(error, 'preserve_workspace_lock', False):
                    lock_state['preserve'] = True
                raise
            finally:
                lock_state['depth'] -= 1
            return

        ownership_token = _acquire_workspace_directory_lock(lock_path)
        depths[lock_key] = {
            'depth': 1,
            'token': ownership_token,
            'preserve': False,
        }
        try:
            yield ownership_token
        except BaseException as error:
            if getattr(error, 'preserve_workspace_lock', False):
                depths[lock_key]['preserve'] = True
            raise
        finally:
            preserve_lock = depths[lock_key]['preserve']
            del depths[lock_key]
            if not preserve_lock:
                _release_workspace_directory_lock(
                    lock_path, ownership_token
                )


def clear_multi_gene_artifacts(intermediate_dir):
    """Remove published Multi-gene artifacts while retaining reusable cache."""
    intermediate_path = Path(intermediate_dir)
    for published_file in (
        intermediate_path / 'multi_genes.fa',
        intermediate_path / 'multi_genes_aligned_oneline.fa',
    ):
        try:
            published_file.unlink()
        except FileNotFoundError:
            pass

    for directory, pattern in (
        (intermediate_path / 'multi_genes', '*.fa'),
        (
            intermediate_path / 'multi_genes_aligned',
            '*_aligned_consensus.fa',
        ),
        (
            intermediate_path / 'multi_genes_aligned_oneline',
            '*_aligned_consensus.fa',
        ),
    ):
        if directory.is_dir():
            for published_file in directory.glob(pattern):
                published_file.unlink()


def clear_single_gene_artifacts(intermediate_dir):
    """Remove a stale Single-gene publication when the category is empty."""
    try:
        (Path(intermediate_dir) / 'single_genes.fa').unlink()
    except FileNotFoundError:
        pass


class MultiGeneAlignmentRunner:
    """Run the multi-gene alignment helper with a content-addressed image."""

    DEFAULT_WORKERS = 4
    DEFAULT_CACHE_ENABLED = True
    DEFAULT_MEMORY = '4g'
    IMAGE_REPOSITORY = 'twheel-run-alignment-and-mask'
    BUILD_CONTEXT = (
        Path(__file__).resolve().parents[1]
        / 'docker'
        / 'run_alignment_and_mask'
    )
    _IGNORED_DIRECTORIES = {
        '.git',
        '.mypy_cache',
        '.pytest_cache',
        '.ruff_cache',
        '__pycache__',
    }
    _IGNORED_FILENAMES = {'.DS_Store'}
    _IGNORED_SUFFIXES = {'.pyc', '.pyo', '.swp', '.tmp'}
    _MEMORY_PATTERN = re.compile(
        r'^(?:[1-9][0-9]*(?:\.[0-9]+)?|0\.[0-9]*[1-9][0-9]*)[bkmg]$',
        re.IGNORECASE,
    )

    def __init__(
        self,
        intermediate_dir,
        config,
        logger=None,
        run_command=subprocess.run,
    ):
        """Initialize paths, command execution, and validated configuration."""
        self.intermediate_dir = Path(intermediate_dir).expanduser().resolve()
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        if not callable(run_command):
            raise TypeError('run_command must be callable')
        self.run_command = run_command

        settings = self._validated_settings(config)
        self.workers = settings['workers']
        self.cache_enabled = settings['cache_enabled']
        self.memory = settings['memory']

        self.build_context = self.BUILD_CONTEXT
        self.input_dir = self.intermediate_dir / 'multi_genes'
        self.output_dir = self.intermediate_dir / 'multi_genes_aligned'
        self.cache_dir = (
            self.intermediate_dir / 'multi_genes_alignment_cache'
        )
        self.buildx_config_dir = Path(tempfile.gettempdir()) / (
            f'twheel-buildx-{os.getuid()}'
        )

    @classmethod
    def _validated_settings(cls, config):
        """Return validated ``multi_gene_alignment`` configuration values."""
        if not isinstance(config, Mapping):
            raise TypeError('config must be a mapping')

        section = config.get('multi_gene_alignment', {})
        if not isinstance(section, Mapping):
            raise TypeError('multi_gene_alignment must be a mapping')

        workers = section.get('workers', cls.DEFAULT_WORKERS)
        if type(workers) is not int or workers <= 0:
            raise ValueError(
                'multi_gene_alignment.workers must be a positive integer'
            )

        cache_enabled = section.get(
            'cache_enabled', cls.DEFAULT_CACHE_ENABLED
        )
        if type(cache_enabled) is not bool:
            raise ValueError(
                'multi_gene_alignment.cache_enabled must be a boolean'
            )

        memory = section.get('memory', cls.DEFAULT_MEMORY)
        if (
            type(memory) is not str
            or not cls._MEMORY_PATTERN.fullmatch(memory)
        ):
            raise ValueError(
                'multi_gene_alignment.memory must be a positive Docker memory '
                'value with a b, k, m, or g suffix'
            )

        return {
            'workers': workers,
            'cache_enabled': cache_enabled,
            'memory': memory,
        }

    def _build_context_digest(self):
        """Hash every relevant build-context path and its contents."""
        context = Path(self.build_context)
        if not context.is_dir():
            raise FileNotFoundError(
                f'Multi-gene alignment build context not found: {context}'
            )

        included_files = []
        for root, directories, filenames in os.walk(context):
            directories[:] = sorted(
                directory
                for directory in directories
                if directory not in self._IGNORED_DIRECTORIES
            )
            root_path = Path(root)
            for filename in sorted(filenames):
                path = root_path / filename
                if not self._include_build_context_file(path):
                    continue
                included_files.append(path)

        if not included_files:
            raise RuntimeError(
                f'Multi-gene alignment build context is empty: {context}'
            )

        digest = hashlib.sha256()
        for path in sorted(
            included_files,
            key=lambda candidate: candidate.relative_to(context).as_posix(),
        ):
            relative_path = path.relative_to(context).as_posix()
            digest.update(relative_path.encode('utf-8'))
            digest.update(b'\0')
            with path.open('rb') as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                    digest.update(chunk)
            digest.update(b'\0')
        return digest.hexdigest()

    def _include_build_context_file(self, path):
        """Return whether a file belongs in the build-context identity."""
        if not path.is_file():
            return False
        if path.name in self._IGNORED_FILENAMES:
            return False
        if path.suffix.lower() in self._IGNORED_SUFFIXES:
            return False
        return not any(
            part in self._IGNORED_DIRECTORIES
            for part in path.parts
        )

    def _inspect_image(self, image_tag):
        """Return the image identity, or ``None`` when the tag is absent."""
        command = [
            'docker',
            'image',
            'inspect',
            '--format',
            '{{.Id}}',
            image_tag,
        ]
        result = self.run_command(
            command,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            return None
        image_identity = (result.stdout or '').strip()
        if not image_identity:
            raise RuntimeError(
                'Docker returned an empty identity while inspecting '
                f'{image_tag}'
            )
        return image_identity

    def _ensure_image(self):
        """Return a built image tag and its immutable tool identity."""
        context_digest = self._build_context_digest()
        image_tag = f'{self.IMAGE_REPOSITORY}:{context_digest}'

        tool_identity = self._inspect_image(image_tag)
        if tool_identity is not None:
            self.logger.info(
                'Reusing multi-gene alignment image: %s', image_tag
            )
            return image_tag, tool_identity

        lock_path = Path('/tmp') / (
            f'twheel-multi-gene-alignment-{context_digest}.lock'
        )
        with _exclusive_file_lock(lock_path):
            tool_identity = self._inspect_image(image_tag)
            if tool_identity is None:
                build_command = [
                    'docker',
                    'build',
                    '-t',
                    image_tag,
                    str(Path(self.build_context)),
                ]
                self.buildx_config_dir.mkdir(
                    parents=True,
                    exist_ok=True,
                )
                build_environment = os.environ.copy()
                # Older releases used sudo and may have left the default
                # buildx state owned by root. Keep tWhEEL build state in a
                # user-writable location without changing Docker context.
                build_environment['BUILDX_CONFIG'] = str(
                    self.buildx_config_dir
                )
                build_result = self.run_command(
                    build_command,
                    capture_output=True,
                    text=True,
                    env=build_environment,
                )
                if build_result.returncode != 0:
                    error = (build_result.stderr or '').strip()
                    raise RuntimeError(
                        'Multi-gene alignment image build failed with '
                        f'exit code {build_result.returncode}: '
                        f'{error or "No error details provided."}'
                    )

                tool_identity = self._inspect_image(image_tag)
                if tool_identity is None:
                    raise RuntimeError(
                        'Built multi-gene alignment image cannot be inspected: '
                        f'{image_tag}'
                    )
                self.logger.info(
                    'Built multi-gene alignment image: %s', image_tag
                )
            else:
                self.logger.info(
                    'Reusing multi-gene alignment image built by another '
                    'process: %s',
                    image_tag,
                )

        return image_tag, tool_identity

    def _normalize_expected_ids(self, expected_multi_ids):
        """Validate and normalize expected Multi-gene identifiers."""
        if expected_multi_ids is None:
            return {
                input_path.stem
                for input_path in self.input_dir.glob('*.fa')
            }
        if isinstance(expected_multi_ids, (str, bytes)):
            raise TypeError('expected_multi_ids must be an iterable of IDs')
        try:
            normalized = set(expected_multi_ids)
        except TypeError as error:
            raise TypeError(
                'expected_multi_ids must be an iterable of IDs'
            ) from error
        invalid_identifiers = (
            type(identifier) is not str or not identifier
            for identifier in normalized
        )
        if any(invalid_identifiers):
            raise ValueError(
                'expected_multi_ids must contain non-empty strings only'
            )
        return normalized

    def _backup_published_outputs(self):
        """Move the current complete output set aside for failure rollback."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        backup_dir = Path(tempfile.mkdtemp(
            dir=self.intermediate_dir,
            prefix='.multi-gene-alignment-backup-',
        ))
        try:
            for published_output in self.output_dir.glob(
                '*_aligned_consensus.fa'
            ):
                backup_path = backup_dir / published_output.name
                os.replace(published_output, backup_path)
        except BaseException as backup_error:
            rollback_errors = []
            for backup_path in backup_dir.glob(
                '*_aligned_consensus.fa'
            ):
                published_output = self.output_dir / backup_path.name
                try:
                    os.replace(backup_path, published_output)
                except BaseException as rollback_error:
                    rollback_errors.append(rollback_error)
            if rollback_errors:
                raise WorkspaceLockSafetyError(
                    'Failed to create the alignment output backup and could '
                    f'not fully restore it. Preserved backup: {backup_dir}'
                ) from backup_error
            shutil.rmtree(backup_dir, ignore_errors=True)
            raise
        return backup_dir

    def _restore_published_outputs(self, backup_dir):
        """Restore the prior complete set after a failed generation."""
        for failed_output in self.output_dir.glob(
            '*_aligned_consensus.fa'
        ):
            failed_output.unlink()
        for previous_output in Path(backup_dir).glob(
            '*_aligned_consensus.fa'
        ):
            os.replace(
                previous_output,
                self.output_dir / previous_output.name,
            )

    def _docker_run_command(
        self,
        image_reference,
        tool_identity,
        parent_lock_token,
        container_name,
    ):
        """Build the Docker command for one alignment run."""
        command = [
            'docker',
            'run',
            '--rm',
            '--name',
            container_name,
            '-m',
            self.memory,
            '-v',
            f'{self.intermediate_dir}:/data',
            '-w',
            '/data',
            image_reference,
            'python',
            '/app/run_alignment_and_mask.py',
            '-i',
            'multi_genes',
            '-o',
            'multi_genes_aligned',
            '--workers',
            str(self.workers),
            '--tool-identity',
            tool_identity,
            '--no-run-lock',
            '--parent-lock-token',
            parent_lock_token,
        ]
        if self.cache_enabled:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            command.extend([
                '--cache-dir',
                'multi_genes_alignment_cache',
            ])
        else:
            command.append('--no-cache')
        return command

    def _stop_failed_container(self, container_name):
        """Stop an interrupted container or preserve the workspace lock."""
        command = ['docker', 'rm', '-f', container_name]
        try:
            result = self.run_command(
                command,
                capture_output=True,
                text=True,
            )
        except BaseException as error:
            raise WorkspaceLockSafetyError(
                'Could not stop or rule out an interrupted Multi-gene '
                f'alignment container: {container_name}. The workspace lock '
                'was preserved; verify Docker state before removing it.'
            ) from error
        error_text = '\n'.join(
            value.strip()
            for value in (result.stdout, result.stderr)
            if value and value.strip()
        )
        if (
            result.returncode != 0
            and 'no such container' not in error_text.lower()
        ):
            raise WorkspaceLockSafetyError(
                'Could not stop or rule out an interrupted Multi-gene '
                f'alignment container: {container_name}. Docker reported: '
                f'{error_text or "no details"}. The workspace lock was '
                'preserved; verify Docker state before removing it.'
            )

    def _validated_output_paths(self, expected_multi_ids):
        """Return outputs after enforcing an exact expected/observed match."""
        output_paths = sorted(
            self.output_dir.glob('*_aligned_consensus.fa')
        )
        suffix = '_aligned_consensus.fa'
        observed_multi_ids = {
            output_path.name[:-len(suffix)]
            for output_path in output_paths
        }
        missing = sorted(expected_multi_ids - observed_multi_ids)
        unexpected = sorted(observed_multi_ids - expected_multi_ids)
        if missing or unexpected:
            details = []
            if missing:
                details.append('missing=' + ','.join(missing))
            if unexpected:
                details.append('unexpected=' + ','.join(unexpected))
            raise RuntimeError(
                'Multi-gene alignment outputs do not match input groups: '
                + '; '.join(details)
            )
        for output_path in output_paths:
            multi_id = output_path.name[:-len(suffix)]
            self._validate_consensus_fasta(output_path, multi_id)
        return [str(path) for path in output_paths]

    @staticmethod
    def _validate_consensus_fasta(output_path, expected_multi_id):
        """Require one non-empty FASTA record with the expected identifier."""
        record_ids = []
        sequence_lengths = []
        current_sequence_length = None
        with Path(output_path).open('r', encoding='utf-8') as handle:
            for line in handle:
                stripped_line = line.strip()
                if not stripped_line:
                    continue
                if stripped_line.startswith('>'):
                    if current_sequence_length is not None:
                        sequence_lengths.append(current_sequence_length)
                    header_text = stripped_line[1:].strip()
                    record_ids.append(
                        header_text.split(maxsplit=1)[0]
                        if header_text else ''
                    )
                    current_sequence_length = 0
                else:
                    if current_sequence_length is None:
                        raise RuntimeError(
                            'Invalid Multi-gene consensus FASTA: '
                            f'{output_path}'
                        )
                    current_sequence_length += len(stripped_line)
        if current_sequence_length is not None:
            sequence_lengths.append(current_sequence_length)

        if (
            record_ids != [expected_multi_id]
            or len(sequence_lengths) != 1
            or sequence_lengths[0] <= 0
        ):
            raise RuntimeError(
                f'Invalid Multi-gene consensus FASTA: {output_path}'
            )

    def run(self, expected_multi_ids):
        """Build the helper if needed, run it, and validate every output."""
        with annotation_workspace_lock(
            self.intermediate_dir
        ) as parent_lock_token:
            if not self.input_dir.is_dir():
                raise FileNotFoundError(
                    'Multi-gene alignment input directory not found: '
                    f'{self.input_dir}'
                )
            expected = self._normalize_expected_ids(expected_multi_ids)
            return self._run_locked(expected, parent_lock_token)

    def _run_locked(self, expected, parent_lock_token):
        """Execute one serialized alignment generation."""
        _, tool_identity = self._ensure_image()
        backup_dir = self._backup_published_outputs()
        remove_backup = False
        validated_paths = None
        container_name = f'twheel-alignment-{uuid.uuid4().hex}'
        try:
            command = self._docker_run_command(
                tool_identity,
                tool_identity,
                parent_lock_token,
                container_name,
            )
            try:
                result = self.run_command(
                    command,
                    capture_output=True,
                    text=True,
                )
            except BaseException:
                self._stop_failed_container(container_name)
                raise
            if result.returncode != 0:
                self._stop_failed_container(container_name)
                error = (result.stderr or '').strip()
                raise RuntimeError(
                    'Multi-gene alignment failed with exit code '
                    f'{result.returncode}: '
                    f'{error or "No error details provided."}'
                )

            summary = (result.stdout or '').strip()
            if summary:
                self.logger.info(
                    'Multi-gene alignment summary:\n%s', summary
                )

            validated_paths = self._validated_output_paths(expected)
        except BaseException as run_error:
            if getattr(
                run_error,
                'preserve_workspace_lock',
                False,
            ):
                raise
            try:
                self._restore_published_outputs(backup_dir)
            except BaseException as rollback_error:
                raise WorkspaceLockSafetyError(
                    'Multi-gene alignment failed and the previous output set '
                    f'could not be fully restored. Preserved backup: {backup_dir}'
                ) from rollback_error
            remove_backup = True
            raise
        else:
            remove_backup = True
        finally:
            if remove_backup:
                shutil.rmtree(backup_dir, ignore_errors=True)
        return validated_paths
