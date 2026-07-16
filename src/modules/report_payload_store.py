"""Immutable browser-payload publication for generated reports."""

import fcntl
import hashlib
import json
import os
import re
import shutil
import tempfile
import threading
import uuid
from collections.abc import Mapping
from contextlib import contextmanager
from pathlib import Path


class ReportPayloadStore:
    """Publish browser payloads as immutable, file-addressable generations."""

    SCHEMA_VERSION = 2
    _NAMESPACE_PATTERN = re.compile(r'^[A-Za-z0-9][A-Za-z0-9_.-]*$')
    _GENERATION_ID_PATTERN = re.compile(r'^[0-9a-f]{32}$')
    _EXTRA_GLOBAL_PATTERN = re.compile(
        r'^__TWHEEL_[A-Za-z0-9_]+__$'
    )
    _RESERVED_GLOBALS = frozenset({
        '__TWHEEL_REPORT_PAYLOAD__',
        '__TWHEEL_REPORT_PAYLOADS__',
        '__TWHEEL_REPORT_GUIDE_CHUNKS__',
    })

    def __init__(self, web_data_dir, guide_chunk_rows=5000):
        """Initialize storage rooted at an output directory's ``web_data``."""
        if (
            isinstance(guide_chunk_rows, bool)
            or not isinstance(guide_chunk_rows, int)
            or guide_chunk_rows <= 0
        ):
            raise ValueError('guide_chunk_rows must be a positive integer')

        self.web_data_dir = Path(web_data_dir).absolute()
        self.root = self.web_data_dir / 'report_payloads'
        self.committed_root = self.root / '.committed'
        self._lock_path = self.root / '.publication.lock'
        self._lock_state = threading.local()
        self.guide_chunk_rows = guide_chunk_rows

    @contextmanager
    def publication_lock(self):
        """Serialize payload and HTML commit operations for one output root."""
        depth = getattr(self._lock_state, 'depth', 0)
        if depth:
            self._lock_state.depth = depth + 1
            try:
                yield
            finally:
                self._lock_state.depth -= 1
            return

        self.root.mkdir(parents=True, exist_ok=True)
        with self._lock_path.open('a+b') as lock_handle:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            self._lock_state.depth = 1
            try:
                yield
            finally:
                self._lock_state.depth = 0
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)

    def publish(
        self,
        namespace,
        core_data,
        guide_rows=(),
        extra_globals=None,
    ):
        """Publish one generation under the output-root publication lock."""
        with self.publication_lock():
            return self._publish_unlocked(
                namespace,
                core_data,
                guide_rows,
                extra_globals,
            )

    def _publish_unlocked(
        self,
        namespace,
        core_data,
        guide_rows,
        extra_globals,
    ):
        """Publish one complete generation and return its browser entry point.

        Guide rows are consumed once and written in bounded chunks. The completed
        temporary directory becomes visible through one atomic directory rename.
        """
        namespace = self._validate_namespace(namespace)
        normalized_globals = self._normalize_extra_globals(extra_globals)
        if guide_rows is None:
            raise TypeError('guide_rows must be an iterable')

        self.root.mkdir(parents=True, exist_ok=True)
        stale_prefix = f'.{namespace}-'
        for candidate in self.root.iterdir():
            if (
                candidate.is_dir()
                and candidate.name.startswith(stale_prefix)
                and '.tmp-' in candidate.name
            ):
                shutil.rmtree(candidate, ignore_errors=True)
        generation_id = uuid.uuid4().hex
        generation_name = f'{namespace}-{generation_id}'
        final_directory = self.root / generation_name
        temporary_directory = Path(tempfile.mkdtemp(
            dir=self.root,
            prefix=f'.{generation_name}.tmp-',
        ))
        guides_directory = temporary_directory / 'guides'

        try:
            guides_directory.mkdir()
            chunks, guide_count, chunk_files = self._write_guide_chunks(
                guides_directory,
                namespace,
                generation_id,
                guide_rows,
            )

            core_payload = {
                'schema_version': self.SCHEMA_VERSION,
                'namespace': namespace,
                'generation_id': generation_id,
                'data': core_data,
                'guides': {
                    'row_count': guide_count,
                    'chunk_rows': self.guide_chunk_rows,
                    'chunks': chunks,
                },
            }
            core_javascript = self._build_core_javascript(
                core_payload,
                normalized_globals,
            )
            core_metadata = self._write_file(
                temporary_directory / 'core.js',
                core_javascript,
            )

            files = {'core.js': core_metadata}
            files.update(chunk_files)
            manifest = {
                'schema_version': self.SCHEMA_VERSION,
                'namespace': namespace,
                'generation_id': generation_id,
                'guide_count': guide_count,
                'guide_chunk_rows': self.guide_chunk_rows,
                'files': files,
            }
            manifest_text = self._json_dumps(manifest, indent=2) + '\n'
            self._write_file(
                temporary_directory / 'manifest.json',
                manifest_text,
            )

            self._fsync_directory(guides_directory)
            self._fsync_directory(temporary_directory)
            os.replace(temporary_directory, final_directory)
            self._fsync_directory(self.root)

            core_path = final_directory / 'core.js'
            manifest_path = final_directory / 'manifest.json'
            core_relative_path = core_path.relative_to(
                self.web_data_dir.parent
            ).as_posix()
            return {
                'generation_id': generation_id,
                'core_path': str(core_path),
                'core_relative_path': core_relative_path,
                'manifest_path': str(manifest_path),
                'guide_count': guide_count,
                'chunks': chunks,
            }
        finally:
            if temporary_directory.exists():
                shutil.rmtree(temporary_directory, ignore_errors=True)

    def cleanup_except(self, namespace, generation_id):
        """Explicitly remove all but one generation under the store lock."""
        namespace = self._validate_namespace(namespace)
        generation_id = self._validate_generation_id(generation_id)
        with self.publication_lock():
            self._cleanup_except_unlocked(namespace, generation_id)

    def _cleanup_except_unlocked(self, namespace, generation_id):
        if not self.root.is_dir():
            return

        keep_name = f'{namespace}-{generation_id}'
        namespace_prefix = f'{namespace}-'
        for candidate in self.root.iterdir():
            if (
                not candidate.is_dir()
                or candidate.name == keep_name
                or not candidate.name.startswith(namespace_prefix)
            ):
                continue

            candidate_generation = candidate.name[len(namespace_prefix):]
            if not self._GENERATION_ID_PATTERN.fullmatch(candidate_generation):
                continue
            try:
                shutil.rmtree(candidate)
            except OSError:
                # Cleanup is best-effort after the new HTML commit is durable.
                continue
            try:
                (self.committed_root / candidate.name).unlink()
            except OSError:
                pass

        self._fsync_directory(self.root)

    def mark_committed(self, namespace, generation_id):
        """Mark a payload after its referring HTML has been committed."""
        namespace = self._validate_namespace(namespace)
        generation_id = self._validate_generation_id(generation_id)
        generation_directory = self.root / f'{namespace}-{generation_id}'
        if not generation_directory.is_dir():
            raise FileNotFoundError(
                f'Report payload generation is missing: {generation_directory}'
            )
        self.committed_root.mkdir(parents=True, exist_ok=True)
        marker_path = self.committed_root / generation_directory.name
        temporary_marker = self.committed_root / (
            f'.{generation_directory.name}.{uuid.uuid4().hex}.tmp'
        )
        try:
            self._write_file(temporary_marker, generation_id + '\n')
            os.replace(temporary_marker, marker_path)
            self._fsync_directory(self.committed_root)
        finally:
            try:
                temporary_marker.unlink()
            except OSError:
                pass

    def cleanup_uncommitted_locked(self, namespace, keep_generation_id):
        """Remove failed publications while the caller holds the store lock."""
        namespace = self._validate_namespace(namespace)
        keep_generation_id = self._validate_generation_id(
            keep_generation_id
        )
        if not self.root.is_dir():
            return

        namespace_prefix = f'{namespace}-'
        keep_name = f'{namespace}-{keep_generation_id}'
        for candidate in self.root.iterdir():
            if (
                not candidate.is_dir()
                or candidate.name == keep_name
                or not candidate.name.startswith(namespace_prefix)
                or (self.committed_root / candidate.name).is_file()
            ):
                continue
            candidate_generation = candidate.name[len(namespace_prefix):]
            if not self._GENERATION_ID_PATTERN.fullmatch(candidate_generation):
                continue
            try:
                shutil.rmtree(candidate)
            except OSError:
                continue
        self._fsync_directory(self.root)

    def _write_guide_chunks(
        self,
        guides_directory,
        namespace,
        generation_id,
        guide_rows,
    ):
        chunks = []
        files = {}
        pending_rows = []
        guide_count = 0

        for row in guide_rows:
            pending_rows.append(row)
            if len(pending_rows) == self.guide_chunk_rows:
                descriptor, metadata = self._write_guide_chunk(
                    guides_directory,
                    namespace,
                    generation_id,
                    len(chunks),
                    pending_rows,
                )
                chunks.append(descriptor)
                files[descriptor['src']] = metadata
                guide_count += len(pending_rows)
                pending_rows = []

        if pending_rows:
            descriptor, metadata = self._write_guide_chunk(
                guides_directory,
                namespace,
                generation_id,
                len(chunks),
                pending_rows,
            )
            chunks.append(descriptor)
            files[descriptor['src']] = metadata
            guide_count += len(pending_rows)

        return chunks, guide_count, files

    def _write_guide_chunk(
        self,
        guides_directory,
        namespace,
        generation_id,
        index,
        rows,
    ):
        file_name = f'{index:06d}.js'
        source = f'guides/{file_name}'
        registry_key = f'{namespace}-{generation_id}:{index:06d}'
        javascript = self._build_guide_javascript(registry_key, rows)
        file_metadata = self._write_file(
            guides_directory / file_name,
            javascript,
        )
        descriptor = {
            'index': index,
            'row_count': len(rows),
            'src': source,
            'registry_key': registry_key,
        }
        manifest_metadata = dict(file_metadata)
        manifest_metadata.update({
            'index': index,
            'row_count': len(rows),
            'registry_key': registry_key,
        })
        return descriptor, manifest_metadata

    def _build_core_javascript(self, payload, extra_globals):
        payload_json = self._json_dumps(payload)
        globals_json = self._json_dumps(extra_globals)
        return (
            '(function () {\n'
            "    'use strict';\n"
            '    var currentScript = document.currentScript;\n'
            '    var baseUrl = currentScript && currentScript.src\n'
            "        ? new URL('.', currentScript.src).href\n"
            "        : '';\n"
            f'    var payload = {payload_json};\n'
            '    payload.base_url = baseUrl;\n'
            '    if (!window.__TWHEEL_REPORT_PAYLOADS__ ||\n'
            '            typeof window.__TWHEEL_REPORT_PAYLOADS__ !== '
            "'object') {\n"
            '        window.__TWHEEL_REPORT_PAYLOADS__ = Object.create(null);\n'
            '    }\n'
            '    window.__TWHEEL_REPORT_PAYLOADS__[payload.generation_id] = '
            'payload;\n'
            '    window.__TWHEEL_REPORT_PAYLOAD__ = payload;\n'
            f'    var extraGlobals = {globals_json};\n'
            '    Object.keys(extraGlobals).forEach(function (name) {\n'
            '        window[name] = extraGlobals[name];\n'
            '    });\n'
            '}());\n'
        )

    def _build_guide_javascript(self, registry_key, rows):
        registry_key_json = self._json_dumps(registry_key)
        rows_json = self._json_dumps(rows)
        return (
            '(function () {\n'
            "    'use strict';\n"
            '    if (!window.__TWHEEL_REPORT_GUIDE_CHUNKS__ ||\n'
            '            typeof window.__TWHEEL_REPORT_GUIDE_CHUNKS__ !== '
            "'object') {\n"
            '        window.__TWHEEL_REPORT_GUIDE_CHUNKS__ = '
            'Object.create(null);\n'
            '    }\n'
            '    window.__TWHEEL_REPORT_GUIDE_CHUNKS__'
            f'[{registry_key_json}] = {rows_json};\n'
            '}());\n'
        )

    @classmethod
    def _validate_namespace(cls, namespace):
        namespace = str(namespace)
        if not cls._NAMESPACE_PATTERN.fullmatch(namespace):
            raise ValueError(f'Invalid report-payload namespace: {namespace!r}')
        return namespace

    @classmethod
    def _validate_generation_id(cls, generation_id):
        generation_id = str(generation_id)
        if not cls._GENERATION_ID_PATTERN.fullmatch(generation_id):
            raise ValueError(
                f'Invalid report-payload generation ID: {generation_id!r}'
            )
        return generation_id

    @classmethod
    def _normalize_extra_globals(cls, extra_globals):
        if extra_globals is None:
            return {}
        if not isinstance(extra_globals, Mapping):
            raise TypeError('extra_globals must be a mapping or None')

        normalized = {}
        for name, value in extra_globals.items():
            if not isinstance(name, str):
                raise ValueError('Extra-global names must be strings')
            if not cls._EXTRA_GLOBAL_PATTERN.fullmatch(name):
                raise ValueError(f'Unsafe report-payload global name: {name!r}')
            if name in cls._RESERVED_GLOBALS:
                raise ValueError(f'Reserved report-payload global name: {name!r}')
            normalized[name] = value
        return normalized

    @staticmethod
    def _json_dumps(value, indent=None):
        serialized = json.dumps(
            value,
            ensure_ascii=False,
            allow_nan=False,
            indent=indent,
            separators=None if indent is not None else (',', ':'),
        )
        return (
            serialized
            .replace('&', r'\u0026')
            .replace('<', r'\u003c')
            .replace('>', r'\u003e')
            .replace('\u2028', r'\u2028')
            .replace('\u2029', r'\u2029')
        )

    @staticmethod
    def _write_file(path, text):
        encoded = text.encode('utf-8')
        with path.open('wb') as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        return {
            'size': len(encoded),
            'sha256': hashlib.sha256(encoded).hexdigest(),
        }

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
