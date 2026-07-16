"""Durable, immutable storage for files needed to regenerate reports."""

import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import uuid
from contextlib import contextmanager
from pathlib import Path, PurePosixPath


class ReportSourceStore:
    """Publish validated report-source snapshots under ``web_data``."""

    MANIFEST_VERSION = 2
    _NAMESPACE_PATTERN = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')

    def __init__(self, web_data_dir):
        self.web_data_dir = Path(web_data_dir)
        self.root = self.web_data_dir / 'report_sources'
        self.snapshots_root = self.root / 'snapshots'
        self.manifest_path = self.root / 'manifest.json'
        self._lock_path = self.root / '.manifest.lock'

    def publish_namespace(self, namespace, sources):
        """Atomically replace one namespace with an immutable snapshot.

        ``sources`` maps logical keys to ``(source_path, relative_path)`` pairs.
        The returned mapping uses logical keys and absolute snapshot paths.
        Publishing an empty mapping removes the namespace from the manifest.
        """
        namespace = self._validate_namespace(namespace)
        normalized_sources = self._normalize_sources(sources)

        if not normalized_sources:
            with self._manifest_lock():
                entries = self._read_manifest_entries()
                entries = self._without_namespace(entries, namespace)
                self._write_manifest_atomic(entries)
                self._cleanup_unreferenced_snapshots(entries)
            return {}

        self.snapshots_root.mkdir(parents=True, exist_ok=True)
        snapshot_name = f'{namespace}-{uuid.uuid4().hex}'
        final_snapshot = self.snapshots_root / snapshot_name
        temporary_snapshot = Path(tempfile.mkdtemp(
            dir=self.snapshots_root,
            prefix=f'.{snapshot_name}.tmp-',
        ))

        namespace_entries = {}
        published_paths = {}
        resolved_published_paths = {}
        snapshot_published = False
        try:
            for logical_key, source, relative_path in normalized_sources:
                destination = temporary_snapshot.joinpath(*relative_path.parts)
                destination.parent.mkdir(parents=True, exist_ok=True)
                self._copy_independent(source, destination)
                size = destination.stat().st_size
                digest = self._hash_file(destination)
                manifest_key = f'{namespace}.{logical_key}'
                final_path = final_snapshot.joinpath(*relative_path.parts)
                published_paths[logical_key] = final_path
                namespace_entries[manifest_key] = {
                    'path': final_path.relative_to(
                        self.web_data_dir
                    ).as_posix(),
                    'size': size,
                    'sha256': digest,
                }

            with self._manifest_lock():
                entries = self._read_manifest_entries()
                updated_entries = self._without_namespace(entries, namespace)
                updated_entries.update(namespace_entries)

                os.replace(temporary_snapshot, final_snapshot)
                snapshot_published = True
                try:
                    self._write_manifest_atomic(updated_entries)
                except Exception:
                    shutil.rmtree(final_snapshot, ignore_errors=True)
                    snapshot_published = False
                    raise

                self._cleanup_unreferenced_snapshots(updated_entries)
                resolved_published_paths = {
                    logical_key: str(path.resolve(strict=True))
                    for logical_key, path in published_paths.items()
                }

            return resolved_published_paths
        finally:
            if temporary_snapshot.exists():
                shutil.rmtree(temporary_snapshot, ignore_errors=True)
            if snapshot_published and not self.manifest_path.exists():
                shutil.rmtree(final_snapshot, ignore_errors=True)

    def load_namespace(self, namespace):
        """Load and validate all published paths for one namespace."""
        namespace = self._validate_namespace(namespace)
        entries = self._read_manifest_entries()
        prefix = f'{namespace}.'
        loaded = {}

        for manifest_key, entry in sorted(entries.items()):
            if not manifest_key.startswith(prefix):
                continue
            logical_key = manifest_key[len(prefix):]
            if not logical_key:
                raise RuntimeError(
                    f'Invalid report-source manifest key: {manifest_key!r}'
                )
            loaded[logical_key] = str(
                self._validate_manifest_entry(namespace, manifest_key, entry)
            )

        return loaded

    @classmethod
    def _validate_namespace(cls, namespace):
        namespace = str(namespace)
        if not cls._NAMESPACE_PATTERN.fullmatch(namespace):
            raise ValueError(f'Invalid report-source namespace: {namespace!r}')
        return namespace

    def _normalize_sources(self, sources):
        if sources is None:
            raise TypeError('sources must be a mapping')

        normalized = []
        destinations = set()
        for logical_key, source_spec in sources.items():
            logical_key = str(logical_key)
            if not logical_key:
                raise ValueError('Report-source logical keys must not be empty')
            try:
                source_path, relative_path = source_spec
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f'Invalid source specification for {logical_key!r}'
                ) from error

            source = Path(source_path)
            if not source.is_file():
                raise FileNotFoundError(
                    f'Report source does not exist: {source}'
                )

            relative = PurePosixPath(str(relative_path))
            if (
                relative.is_absolute()
                or not relative.parts
                or '..' in relative.parts
            ):
                raise ValueError(
                    f'Invalid report-source relative path: {relative_path!r}'
                )
            if relative in destinations:
                raise ValueError(
                    f'Duplicate report-source destination: {relative}'
                )
            destinations.add(relative)
            normalized.append((logical_key, source, relative))

        return normalized

    def _validate_manifest_entry(self, namespace, manifest_key, entry):
        if not isinstance(entry, dict):
            raise RuntimeError(
                f'Invalid report-source manifest entry: {manifest_key!r}'
            )

        relative_path = entry.get('path')
        expected_size = entry.get('size')
        expected_hash = entry.get('sha256')
        if (
            not isinstance(relative_path, str)
            or not isinstance(expected_size, int)
            or expected_size < 0
            or not isinstance(expected_hash, str)
            or not re.fullmatch(r'[0-9a-f]{64}', expected_hash)
        ):
            raise RuntimeError(
                f'Invalid report-source manifest metadata: {manifest_key!r}'
            )

        relative = PurePosixPath(relative_path)
        if relative.is_absolute() or not relative.parts or '..' in relative.parts:
            raise RuntimeError(
                f'Report-source path escapes web_data: {relative_path!r}'
            )

        candidate = self.web_data_dir.joinpath(*relative.parts)
        try:
            resolved_candidate = candidate.resolve(strict=True)
        except OSError as error:
            raise RuntimeError(
                f'Report source is missing: {relative_path!r}'
            ) from error

        resolved_snapshots_root = self.snapshots_root.resolve(strict=False)
        try:
            common_path = os.path.commonpath([
                str(resolved_candidate),
                str(resolved_snapshots_root),
            ])
        except ValueError as error:
            raise RuntimeError(
                f'Report-source path escapes snapshot storage: {relative_path!r}'
            ) from error
        if common_path != str(resolved_snapshots_root):
            raise RuntimeError(
                f'Report-source path escapes snapshot storage: {relative_path!r}'
            )

        try:
            snapshot_relative = resolved_candidate.relative_to(
                resolved_snapshots_root
            )
        except ValueError as error:
            raise RuntimeError(
                f'Report-source path escapes snapshot storage: {relative_path!r}'
            ) from error
        if (
            len(snapshot_relative.parts) < 2
            or not snapshot_relative.parts[0].startswith(f'{namespace}-')
        ):
            raise RuntimeError(
                f'Report source is outside namespace {namespace!r}: '
                f'{relative_path!r}'
            )
        if not resolved_candidate.is_file():
            raise RuntimeError(f'Report source is not a file: {relative_path!r}')

        actual_size = resolved_candidate.stat().st_size
        if actual_size != expected_size:
            raise RuntimeError(
                f'Report-source size mismatch for {manifest_key!r}: '
                f'expected {expected_size}, found {actual_size}'
            )
        actual_hash = self._hash_file(resolved_candidate)
        if actual_hash != expected_hash:
            raise RuntimeError(
                f'Report-source checksum mismatch for {manifest_key!r}'
            )
        return resolved_candidate

    def _read_manifest_entries(self):
        if not self.manifest_path.is_file():
            return {}
        try:
            with self.manifest_path.open(encoding='utf-8') as handle:
                payload = json.load(handle)
        except (OSError, ValueError, TypeError) as error:
            raise RuntimeError('Cannot read the report-source manifest') from error

        if payload.get('version') != self.MANIFEST_VERSION:
            raise RuntimeError(
                'Unsupported report-source manifest version: '
                f'{payload.get("version")!r}'
            )
        entries = payload.get('sources')
        if not isinstance(entries, dict):
            raise RuntimeError('Invalid report-source manifest sources')
        return dict(entries)

    def _write_manifest_atomic(self, entries):
        self.root.mkdir(parents=True, exist_ok=True)
        payload = {
            'version': self.MANIFEST_VERSION,
            'sources': dict(sorted(entries.items())),
        }
        descriptor, temporary_path = tempfile.mkstemp(
            dir=self.root,
            prefix='.manifest.',
            suffix='.tmp',
        )
        try:
            with os.fdopen(descriptor, 'w', encoding='utf-8') as handle:
                json.dump(payload, handle, ensure_ascii=False, indent=2)
                handle.write('\n')
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary_path, self.manifest_path)
            self._fsync_directory(self.root)
        finally:
            if os.path.exists(temporary_path):
                os.unlink(temporary_path)

    @contextmanager
    def _manifest_lock(self):
        self.root.mkdir(parents=True, exist_ok=True)
        with self._lock_path.open('a+b') as lock_handle:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)

    def _cleanup_unreferenced_snapshots(self, entries):
        if not self.snapshots_root.is_dir():
            return

        referenced_names = set()
        for entry in entries.values():
            if not isinstance(entry, dict) or not isinstance(entry.get('path'), str):
                continue
            path = PurePosixPath(entry['path'])
            expected_prefix = PurePosixPath(
                self.root.relative_to(self.web_data_dir).as_posix()
            ) / 'snapshots'
            try:
                snapshot_relative = path.relative_to(expected_prefix)
            except ValueError:
                continue
            if snapshot_relative.parts:
                referenced_names.add(snapshot_relative.parts[0])

        for snapshot_path in self.snapshots_root.iterdir():
            if (
                not snapshot_path.is_dir()
                or snapshot_path.name.startswith('.')
                or snapshot_path.name in referenced_names
            ):
                continue
            try:
                shutil.rmtree(snapshot_path)
            except OSError:
                pass

    @staticmethod
    def _without_namespace(entries, namespace):
        prefix = f'{namespace}.'
        return {
            key: value for key, value in entries.items()
            if not key.startswith(prefix)
        }

    @classmethod
    def _copy_independent(cls, source, destination):
        copied = False
        command = None
        if sys.platform == 'darwin':
            command = ['cp', '-c', '-p', str(source), str(destination)]
        elif sys.platform.startswith('linux'):
            command = [
                'cp', '--reflink=always', '--preserve=mode,timestamps',
                str(source), str(destination),
            ]

        if command is not None:
            try:
                result = subprocess.run(
                    command,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    check=False,
                )
                copied = result.returncode == 0 and destination.is_file()
            except OSError:
                copied = False

        if copied:
            try:
                if os.path.samefile(source, destination):
                    copied = False
            except OSError:
                copied = False

        if not copied:
            if destination.exists() or destination.is_symlink():
                destination.unlink()
            shutil.copy2(source, destination)

        if os.path.samefile(source, destination):
            destination.unlink()
            raise RuntimeError(
                f'Report source snapshot shares an inode with its input: {source}'
            )

    @staticmethod
    def _hash_file(file_path):
        digest = hashlib.sha256()
        with open(file_path, 'rb') as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(chunk)
        return digest.hexdigest()

    @staticmethod
    def _fsync_directory(directory):
        try:
            descriptor = os.open(directory, os.O_RDONLY)
        except OSError:
            return
        try:
            os.fsync(descriptor)
        except OSError:
            pass
        finally:
            os.close(descriptor)
