"""Validated, module-specific cache handling for FASTA search results."""

import json
import os
import subprocess
import tempfile


class FastaSearchCache:
    """Run a FASTA search and reuse only a matching completed result."""

    CACHE_VERSION = 1

    def __init__(self, intermediate_dir, namespace, logger):
        self.intermediate_dir = os.path.abspath(intermediate_dir)
        self.namespace = namespace
        self.logger = logger
        self.result_path = os.path.join(
            self.intermediate_dir,
            f'ggsearch_results_{namespace}.txt',
        )
        self.metadata_path = os.path.join(
            self.intermediate_dir,
            f'ggsearch_results_{namespace}.meta.json',
        )

    @staticmethod
    def _file_signature(file_path):
        """Return a lightweight signature without hashing large FASTA files."""
        absolute_path = os.path.realpath(file_path)
        stat_result = os.stat(absolute_path)
        return {
            'path': absolute_path,
            'size': stat_result.st_size,
            'mtime_ns': stat_result.st_mtime_ns,
        }

    def _build_metadata(self, query_path, db_path, search_settings):
        """Build the complete cache identity for one FASTA search."""
        return {
            'cache_version': self.CACHE_VERSION,
            'namespace': self.namespace,
            'query': self._file_signature(query_path),
            'database': self._file_signature(db_path),
            'search_settings': search_settings,
        }

    def _is_valid(self, expected_metadata):
        """Return whether result and metadata form a complete matching cache."""
        if not os.path.isfile(self.result_path):
            return False
        if os.path.getsize(self.result_path) == 0:
            return False
        if not os.path.isfile(self.metadata_path):
            return False

        try:
            with open(self.metadata_path, 'r', encoding='utf-8') as metadata_file:
                cached_metadata = json.load(metadata_file)
        except (OSError, json.JSONDecodeError) as error:
            self.logger.warning(
                f"Unable to read FASTA search cache metadata "
                f"{self.metadata_path}: {error}"
            )
            return False

        return cached_metadata == expected_metadata

    def _write_metadata(self, metadata):
        """Atomically replace cache metadata after a successful search."""
        metadata_temp_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode='w',
                encoding='utf-8',
                dir=self.intermediate_dir,
                prefix=f'.ggsearch_results_{self.namespace}.',
                suffix='.meta.tmp',
                delete=False,
            ) as metadata_file:
                metadata_temp_path = metadata_file.name
                json.dump(metadata, metadata_file, indent=2, sort_keys=True)
                metadata_file.write('\n')

            os.replace(metadata_temp_path, self.metadata_path)
            metadata_temp_path = None
        finally:
            if metadata_temp_path and os.path.exists(metadata_temp_path):
                os.unlink(metadata_temp_path)

    def run(self, cmd, query_path, db_path, search_settings):
        """Return a valid cached result or execute the search atomically."""
        expected_metadata = self._build_metadata(
            query_path,
            db_path,
            search_settings,
        )

        if self._is_valid(expected_metadata):
            self.logger.info(
                f"Reusing validated FASTA search cache: {self.result_path}"
            )
            return self.result_path

        if os.path.exists(self.result_path) or os.path.exists(self.metadata_path):
            self.logger.warning(
                f"FASTA search cache is stale or incomplete; rerunning: "
                f"{self.result_path}"
            )

        result_temp_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode='w',
                encoding='utf-8',
                dir=self.intermediate_dir,
                prefix=f'.ggsearch_results_{self.namespace}.',
                suffix='.result.tmp',
                delete=False,
            ) as result_file:
                result_temp_path = result_file.name
                process_result = subprocess.run(
                    cmd,
                    stdout=result_file,
                    stderr=subprocess.PIPE,
                    text=True,
                )

            if process_result.returncode != 0:
                error_message = (
                    (process_result.stderr or '').strip()
                    or 'No error details provided.'
                )
                raise RuntimeError(
                    f"FASTA search failed with exit code "
                    f"{process_result.returncode}: {error_message}"
                )

            if os.path.getsize(result_temp_path) == 0:
                raise RuntimeError("FASTA search produced an empty result.")

            os.replace(result_temp_path, self.result_path)
            result_temp_path = None
            self._write_metadata(expected_metadata)
        finally:
            if result_temp_path and os.path.exists(result_temp_path):
                os.unlink(result_temp_path)

        self.logger.info(f"Saved FASTA search results: {self.result_path}")
        return self.result_path
