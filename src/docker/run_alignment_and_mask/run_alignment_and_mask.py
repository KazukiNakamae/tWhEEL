#!/usr/bin/env python3
"""Align Multi-gene transcript sequences and generate consensus sequences."""

import argparse
import fcntl
import hashlib
import json
import os
import shutil
import socket
import struct
import subprocess
import sys
import tempfile
import time
import uuid
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass(frozen=True)
class GroupDescriptor:
    """Describe source FASTA files for one Multi-gene group."""

    source_paths: tuple
    legacy_sequences: tuple = ()
    strict_single_file: bool = False


class WorkspaceLockSafetyError(RuntimeError):
    """Signal that a standalone workspace lock must remain fail-closed."""

    preserve_workspace_lock = True


class AlignmentProcessor:
    """Align transcript groups and generate one consensus per Multi gene."""

    ALGORITHM_SCHEMA = "majority-nongap-consensus-v1"
    CACHE_FORMAT_VERSION = 1
    DISCOVERY_PROGRESS_INTERVAL = 1000
    WORKSPACE_LOCK_TIMEOUT_SECONDS = 30.0

    def __init__(
        self,
        input_dir,
        output_dir,
        cache_dir=None,
        workers=1,
        tool_identity=None,
        cache_enabled=True,
        run_lock_path=None,
        acquire_run_lock=True,
        parent_lock_token=None,
    ):
        """Initialize the processor with optional caching and parallelism.

        Args:
            input_dir: Directory containing grouped FASTA inputs.
            output_dir: Directory receiving materialized consensus FASTAs.
            cache_dir: Persistent group-cache directory. When omitted, a sibling
                directory next to ``output_dir`` is used.
            workers: Maximum number of groups processed concurrently.
            tool_identity: Caller-supplied tool or image identity included in
                every cache key.
            cache_enabled: Disable only when fresh computation is required.
            run_lock_path: Cross-process workspace lock shared with callers.
            acquire_run_lock: Disable only when the parent already holds it.
            parent_lock_token: Ownership token required when locking is skipped.
        """
        self.input_dir = Path(input_dir).expanduser().resolve()
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if isinstance(workers, bool):
            raise ValueError("workers must be a positive integer")
        try:
            parsed_workers = int(workers)
        except (TypeError, ValueError) as error:
            raise ValueError("workers must be a positive integer") from error
        if parsed_workers < 1:
            raise ValueError("workers must be a positive integer")
        self.workers = parsed_workers

        self.cache_enabled = bool(cache_enabled)
        self.acquire_run_lock = bool(acquire_run_lock)
        self.parent_lock_token = str(parent_lock_token or "").strip()
        if not self.acquire_run_lock and not self.parent_lock_token:
            raise ValueError(
                "parent_lock_token is required when acquire_run_lock is false"
            )
        self.run_lock_path = (
            Path(run_lock_path)
            if run_lock_path is not None
            else self.output_dir.parent / ".annotation-workspace.lock.d"
        ).expanduser().resolve()
        if tool_identity is None:
            tool_identity = (
                self._automatic_tool_identity()
                if self.cache_enabled
                else "cache-disabled"
            )
        self.tool_identity = str(tool_identity or "").strip()
        if not self.tool_identity:
            raise ValueError("tool_identity must be non-empty")

        self.cache_dir = (
            Path(cache_dir).expanduser().resolve()
            if cache_dir is not None
            else self.output_dir.parent / ".run_alignment_and_mask_cache"
        )
        self.cache_objects_dir = self.cache_dir / "objects"
        self.cache_metadata_dir = self.cache_dir / "metadata"
        self.cache_locks_dir = self.cache_dir / "locks"
        if self.cache_enabled:
            for directory in (
                self.cache_objects_dir,
                self.cache_metadata_dir,
                self.cache_locks_dir,
            ):
                directory.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def _automatic_tool_identity():
        """Derive a cache identity from source and installed tool versions."""
        try:
            version_result = subprocess.run(
                ["mafft", "--version"],
                capture_output=True,
                text=True,
            )
        except OSError as error:
            raise RuntimeError(
                "Cannot derive cache identity because MAFFT is unavailable."
            ) from error
        mafft_version = "\n".join(
            part.strip()
            for part in (version_result.stdout, version_result.stderr)
            if part and part.strip()
        )
        if version_result.returncode != 0 or not mafft_version:
            raise RuntimeError(
                "Cannot derive cache identity from `mafft --version`."
            )

        digest = hashlib.sha256()
        for value in (
            Path(__file__).read_bytes(),
            sys.version.encode("utf-8"),
            Bio.__version__.encode("utf-8"),
            mafft_version.encode("utf-8"),
        ):
            digest.update(struct.pack(">Q", len(value)))
            digest.update(value)
        return f"auto-sha256:{digest.hexdigest()}"

    @staticmethod
    def _multi_id_from_record_id(record_id):
        """Return a Multi group ID encoded in a FASTA record identifier."""
        header_parts = str(record_id).split("_")
        if len(header_parts) >= 2 and header_parts[0] == "Multi":
            return f"{header_parts[0]}_{header_parts[1]}"
        return None

    def group_sequences_by_multi_id(self):
        """Discover groups without retaining every sequence in memory."""
        group_sources = {}
        legacy_group_sequences = {}
        strict_groups = set()
        fasta_files = sorted(self.input_dir.glob("*.fa"))
        total_files = len(fasta_files)

        for file_index, fasta_file in enumerate(fasta_files, 1):
            stem_parts = fasta_file.stem.split("_")
            standard_multi_id = (
                fasta_file.stem
                if (
                    len(stem_parts) == 2
                    and stem_parts[0] == "Multi"
                    and bool(stem_parts[1])
                )
                else None
            )
            if standard_multi_id:
                group_sources.setdefault(standard_multi_id, []).append(
                    fasta_file
                )
                strict_groups.add(standard_multi_id)
            else:
                with fasta_file.open("r", encoding="utf-8") as fasta_handle:
                    for record in SeqIO.parse(fasta_handle, "fasta"):
                        multi_id = self._multi_id_from_record_id(record.id)
                        if multi_id:
                            legacy_group_sequences.setdefault(
                                multi_id, []
                            ).append(record)
                            strict_groups.discard(multi_id)

            if (
                file_index % self.DISCOVERY_PROGRESS_INTERVAL == 0
                or file_index == total_files
            ):
                print(
                    "Group discovery progress: "
                    f"{file_index}/{total_files} FASTA files, "
                    f"{len(set(group_sources) | set(legacy_group_sequences))} "
                    "groups",
                    flush=True,
                )

        all_multi_ids = set(group_sources) | set(legacy_group_sequences)
        return {
            multi_id: GroupDescriptor(
                tuple(group_sources.get(multi_id, ())),
                tuple(legacy_group_sequences.get(multi_id, ())),
                strict_single_file=(
                    multi_id in strict_groups
                    and len(group_sources.get(multi_id, ())) == 1
                    and multi_id not in legacy_group_sequences
                ),
            )
            for multi_id in sorted(all_multi_ids)
        }

    def _load_group_sequences(self, multi_id, group_source):
        """Load one group from a descriptor or a legacy sequence collection."""
        if isinstance(group_source, GroupDescriptor):
            sequences = list(group_source.legacy_sequences)
            for source_path in group_source.source_paths:
                with source_path.open("r", encoding="utf-8") as fasta_handle:
                    for record in SeqIO.parse(fasta_handle, "fasta"):
                        record_multi_id = self._multi_id_from_record_id(
                            record.id
                        )
                        if group_source.strict_single_file:
                            if record_multi_id != multi_id:
                                raise RuntimeError(
                                    "Standard group FASTA contains a record "
                                    f"for {record_multi_id or 'an invalid group'}; "
                                    f"expected {multi_id}: {source_path}"
                                )
                            sequences.append(record)
                        elif record_multi_id == multi_id:
                            sequences.append(record)
        else:
            sequences = list(group_source)

        if not sequences:
            raise RuntimeError(f"No sequences were found for {multi_id}.")
        empty_record_ids = [
            record.id for record in sequences if len(record.seq) == 0
        ]
        if empty_record_ids:
            raise RuntimeError(
                f"Empty input sequence(s) for {multi_id}: "
                + ", ".join(empty_record_ids[:5])
            )
        return sequences

    def run_mafft_alignment(self, sequences, multi_id):
        """Run MAFFT and return a validated multiple-sequence alignment."""
        if len(sequences) < 2:
            return sequences

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fa", delete=False
        ) as temp_input:
            temp_input_path = temp_input.name
            SeqIO.write(sequences, temp_input, "fasta")

        temp_output_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fa", delete=False
            ) as temp_output:
                temp_output_path = temp_output.name

            cmd = [
                "mafft",
                "--auto",
                "--quiet",
                "--thread",
                "1",
                temp_input_path,
            ]
            with open(temp_output_path, "w") as output_file:
                result = subprocess.run(
                    cmd,
                    stdout=output_file,
                    stderr=subprocess.PIPE,
                    text=True,
                )

            if result.returncode != 0:
                error_message = (
                    (result.stderr or "").strip()
                    or "No error details provided."
                )
                raise RuntimeError(
                    f"MAFFT failed for {multi_id} with exit code "
                    f"{result.returncode}: {error_message}"
                )

            with open(temp_output_path, "r") as alignment_file:
                aligned_sequences = list(
                    SeqIO.parse(alignment_file, "fasta")
                )
            if not aligned_sequences:
                raise RuntimeError(
                    f"MAFFT produced no aligned sequences for {multi_id}."
                )
            if len(aligned_sequences) != len(sequences):
                raise RuntimeError(
                    f"MAFFT returned {len(aligned_sequences)} of "
                    f"{len(sequences)} sequences for {multi_id}."
                )

            alignment_lengths = {
                len(record.seq) for record in aligned_sequences
            }
            if len(alignment_lengths) != 1 or 0 in alignment_lengths:
                raise RuntimeError(
                    f"MAFFT produced an invalid alignment for {multi_id}: "
                    f"sequence lengths={sorted(alignment_lengths)}."
                )
            return aligned_sequences
        finally:
            if os.path.exists(temp_input_path):
                os.unlink(temp_input_path)
            if temp_output_path and os.path.exists(temp_output_path):
                os.unlink(temp_output_path)

    def generate_consensus(self, aligned_sequences, multi_id):
        """Generate a consensus sequence from a validated alignment."""
        if not aligned_sequences:
            return None

        if len(aligned_sequences) == 1:
            return SeqRecord(
                aligned_sequences[0].seq,
                id=multi_id,
                description="consensus sequence",
            )

        alignment_lengths = {
            len(record.seq) for record in aligned_sequences
        }
        if len(alignment_lengths) != 1 or 0 in alignment_lengths:
            raise ValueError(
                f"Consensus input is not a valid alignment for {multi_id}: "
                f"sequence lengths={sorted(alignment_lengths)}."
            )

        alignment_length = len(aligned_sequences[0])
        consensus_sequence = []
        for position in range(alignment_length):
            character_counts = {}
            for sequence_record in aligned_sequences:
                character = str(sequence_record.seq[position]).upper()
                if character != "-":
                    character_counts[character] = (
                        character_counts.get(character, 0) + 1
                    )

            if character_counts:
                consensus_sequence.append(
                    max(character_counts, key=character_counts.get)
                )
            else:
                consensus_sequence.append("N")

        consensus_seq = Seq("".join(consensus_sequence).replace("-", ""))
        return SeqRecord(
            consensus_seq,
            id=multi_id,
            description=(
                "consensus sequence from "
                f"{len(aligned_sequences)} sequences"
            ),
        )

    def _compute_consensus_record(self, multi_id, sequences):
        """Compute one consensus without publishing it."""
        aligned_sequences = self.run_mafft_alignment(sequences, multi_id)
        consensus_record = self.generate_consensus(
            aligned_sequences, multi_id
        )
        if consensus_record is None:
            raise RuntimeError(
                f"Consensus generation returned no result for {multi_id}."
            )
        if consensus_record.id != multi_id or len(consensus_record.seq) == 0:
            raise RuntimeError(
                f"Consensus generation returned an invalid result for {multi_id}."
            )
        return consensus_record

    @staticmethod
    def _fsync_directory(directory):
        """Best-effort directory sync after an atomic replacement."""
        try:
            directory_fd = os.open(directory, os.O_RDONLY)
        except OSError:
            return
        try:
            os.fsync(directory_fd)
        except OSError:
            pass
        finally:
            os.close(directory_fd)

    @staticmethod
    def _write_fasta_temp(record, directory, prefix):
        """Write and sync one FASTA record to a temporary file."""
        directory.mkdir(parents=True, exist_ok=True)
        temporary_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                newline="\n",
                dir=directory,
                prefix=prefix,
                suffix=".tmp",
                delete=False,
            ) as temporary_file:
                temporary_path = Path(temporary_file.name)
                SeqIO.write(record, temporary_file, "fasta")
                temporary_file.flush()
                os.fsync(temporary_file.fileno())
            return temporary_path
        except Exception:
            if temporary_path and temporary_path.exists():
                temporary_path.unlink()
            raise

    def _atomic_write_consensus(self, record, output_path):
        """Publish a final consensus FASTA atomically."""
        output_path = Path(output_path)
        temporary_path = self._write_fasta_temp(
            record,
            output_path.parent,
            f".{output_path.name}.",
        )
        try:
            os.replace(temporary_path, output_path)
            self._fsync_directory(output_path.parent)
        finally:
            if temporary_path.exists():
                temporary_path.unlink()
        return output_path

    def process_multi_gene_group(self, multi_id, sequences):
        """Align one Multi-gene group and atomically save its consensus."""
        sequences = self._load_group_sequences(multi_id, sequences)
        consensus_record = self._compute_consensus_record(multi_id, sequences)
        output_path = (
            self.output_dir / f"{multi_id}_aligned_consensus.fa"
        )
        return self._atomic_write_consensus(consensus_record, output_path)

    @staticmethod
    def _update_framed_hash(digest, value):
        """Add one unambiguous UTF-8 value to a streaming digest."""
        encoded = str(value).encode("utf-8")
        digest.update(struct.pack(">Q", len(encoded)))
        digest.update(encoded)

    def _group_cache_key(self, multi_id, sequences):
        """Hash ordered IDs and sequences with algorithm and tool identity."""
        digest = hashlib.sha256()
        for value in (
            self.CACHE_FORMAT_VERSION,
            self.ALGORITHM_SCHEMA,
            self.tool_identity,
            multi_id,
            len(sequences),
        ):
            self._update_framed_hash(digest, value)
        for record in sequences:
            self._update_framed_hash(digest, record.id)
            self._update_framed_hash(digest, str(record.seq))
        return digest.hexdigest()

    def _cache_paths(self, cache_key):
        """Return object, metadata, and lock paths for one cache key."""
        return (
            self.cache_objects_dir / f"{cache_key}.fa",
            self.cache_metadata_dir / f"{cache_key}.json",
            self.cache_locks_dir / f"{cache_key}.lock",
        )

    @contextmanager
    def _exclusive_file_lock(self, lock_path):
        """Serialize a filesystem publication scope across processes."""
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        with open(lock_path, "a+b") as lock_file:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)

    @contextmanager
    def _workspace_directory_lock(self):
        """Acquire a bind-mount-safe workspace lock with atomic mkdir."""
        lock_dir = Path(self.run_lock_path)
        lock_dir.parent.mkdir(parents=True, exist_ok=True)
        deadline = time.monotonic() + self.WORKSPACE_LOCK_TIMEOUT_SECONDS
        while True:
            try:
                lock_dir.mkdir()
                break
            except FileExistsError:
                if time.monotonic() >= deadline:
                    owner_path = lock_dir / "owner.json"
                    try:
                        owner_details = owner_path.read_text(
                            encoding="utf-8"
                        ).strip()
                    except OSError:
                        owner_details = "owner metadata unavailable"
                    raise RuntimeError(
                        "Timed out waiting for the tWhEEL workspace lock. "
                        f"Lock directory: {lock_dir}. Owner: {owner_details}. "
                        "If the owning process or container crashed, verify "
                        "that no tWhEEL container is still running before "
                        "removing the lock directory manually."
                    )
                time.sleep(0.1)

        ownership_token = uuid.uuid4().hex
        owner_path = lock_dir / "owner.json"
        try:
            self._atomic_write_json(
                {
                    "token": ownership_token,
                    "pid": os.getpid(),
                    "hostname": socket.gethostname(),
                    "created_at_epoch": time.time(),
                },
                owner_path,
            )
        except BaseException:
            shutil.rmtree(lock_dir, ignore_errors=True)
            raise

        release_lock = True
        try:
            yield
        except BaseException as error:
            if getattr(error, "preserve_workspace_lock", False):
                release_lock = False
            raise
        finally:
            if release_lock:
                try:
                    owner_payload = json.loads(
                        owner_path.read_text(encoding="utf-8")
                    )
                except (OSError, ValueError, TypeError) as error:
                    raise RuntimeError(
                        "Cannot verify ownership of workspace lock: "
                        f"{lock_dir}"
                    ) from error
                if owner_payload.get("token") != ownership_token:
                    raise RuntimeError(
                        "Workspace lock ownership changed unexpectedly: "
                        f"{lock_dir}"
                    )
                owner_path.unlink()
                lock_dir.rmdir()

    def _validate_parent_workspace_lock(self):
        """Verify the parent still owns the cross-kernel workspace lock."""
        owner_path = Path(self.run_lock_path) / "owner.json"
        try:
            owner_payload = json.loads(
                owner_path.read_text(encoding="utf-8")
            )
        except (OSError, ValueError, TypeError) as error:
            raise RuntimeError(
                "Parent workspace lock is missing or invalid: "
                f"{self.run_lock_path}"
            ) from error
        if owner_payload.get("token") != self.parent_lock_token:
            raise RuntimeError(
                "Parent workspace lock ownership does not match this "
                "alignment generation."
            )

    @staticmethod
    def _file_sha256(file_path):
        """Return a streaming SHA256 digest for a file."""
        digest = hashlib.sha256()
        with open(file_path, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def _load_valid_cache(self, cache_key, multi_id):
        """Return a validated cache object or ``None`` on any inconsistency."""
        object_path, metadata_path, _ = self._cache_paths(cache_key)
        if not object_path.is_file() or not metadata_path.is_file():
            return None

        try:
            with open(metadata_path, "r", encoding="utf-8") as handle:
                metadata = json.load(handle)
            expected_values = {
                "cache_format_version": self.CACHE_FORMAT_VERSION,
                "cache_key": cache_key,
                "algorithm_schema": self.ALGORITHM_SCHEMA,
                "tool_identity": self.tool_identity,
                "record_id": multi_id,
            }
            if any(
                metadata.get(field) != value
                for field, value in expected_values.items()
            ):
                return None

            output_size = metadata.get("output_size")
            output_sha256 = metadata.get("output_sha256")
            sequence_length = metadata.get("sequence_length")
            if (
                not isinstance(output_size, int)
                or output_size <= 0
                or not isinstance(output_sha256, str)
                or len(output_sha256) != 64
                or not isinstance(sequence_length, int)
                or sequence_length <= 0
            ):
                return None
            if object_path.stat().st_size != output_size:
                return None
            if self._file_sha256(object_path) != output_sha256:
                return None

            with object_path.open("r", encoding="utf-8") as cache_handle:
                cached_records = list(SeqIO.parse(cache_handle, "fasta"))
            if (
                len(cached_records) != 1
                or cached_records[0].id != multi_id
                or len(cached_records[0].seq) != sequence_length
            ):
                return None
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            return None
        return object_path

    def _discard_cache_entry(self, cache_key):
        """Remove incomplete or corrupt cache publications for one key."""
        object_path, metadata_path, _ = self._cache_paths(cache_key)
        for cache_path in (metadata_path, object_path):
            try:
                cache_path.unlink()
            except FileNotFoundError:
                pass

    def _atomic_write_json(self, payload, destination):
        """Publish JSON metadata atomically."""
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="w",
                encoding="utf-8",
                dir=destination.parent,
                prefix=f".{destination.name}.",
                suffix=".tmp",
                delete=False,
            ) as temporary_file:
                temporary_path = Path(temporary_file.name)
                json.dump(payload, temporary_file, sort_keys=True, indent=2)
                temporary_file.write("\n")
                temporary_file.flush()
                os.fsync(temporary_file.fileno())
            os.replace(temporary_path, destination)
            temporary_path = None
            self._fsync_directory(destination.parent)
        finally:
            if temporary_path and temporary_path.exists():
                temporary_path.unlink()

    def _publish_cache_record(self, cache_key, multi_id, consensus_record):
        """Atomically publish a cache object followed by validating metadata."""
        object_path, metadata_path, _ = self._cache_paths(cache_key)
        temporary_object = self._write_fasta_temp(
            consensus_record,
            object_path.parent,
            f".{object_path.name}.",
        )
        try:
            output_size = temporary_object.stat().st_size
            output_sha256 = self._file_sha256(temporary_object)
            os.replace(temporary_object, object_path)
            self._fsync_directory(object_path.parent)
        finally:
            if temporary_object.exists():
                temporary_object.unlink()

        metadata = {
            "cache_format_version": self.CACHE_FORMAT_VERSION,
            "cache_key": cache_key,
            "algorithm_schema": self.ALGORITHM_SCHEMA,
            "tool_identity": self.tool_identity,
            "record_id": multi_id,
            "sequence_length": len(consensus_record.seq),
            "output_size": output_size,
            "output_sha256": output_sha256,
        }
        self._atomic_write_json(metadata, metadata_path)
        return object_path

    def _atomic_materialize_cache(self, cache_object, output_path):
        """Atomically copy a cache object to an independent final file."""
        cache_object = Path(cache_object)
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path = None
        try:
            with tempfile.NamedTemporaryFile(
                mode="wb",
                dir=output_path.parent,
                prefix=f".{output_path.name}.",
                suffix=".tmp",
                delete=False,
            ) as temporary_file:
                temporary_path = Path(temporary_file.name)
                with open(cache_object, "rb") as source_file:
                    shutil.copyfileobj(
                        source_file, temporary_file, length=1024 * 1024
                    )
                temporary_file.flush()
                os.fsync(temporary_file.fileno())
            os.replace(temporary_path, output_path)
            temporary_path = None
            self._fsync_directory(output_path.parent)
        finally:
            if temporary_path and temporary_path.exists():
                temporary_path.unlink()
        return output_path

    def _process_group(self, multi_id, group_source):
        """Compute or reuse one group, then materialize its final output."""
        sequences = self._load_group_sequences(multi_id, group_source)
        output_path = (
            self.output_dir / f"{multi_id}_aligned_consensus.fa"
        )

        if not self.cache_enabled:
            return self.process_multi_gene_group(
                multi_id, sequences
            ), False

        cache_key = self._group_cache_key(multi_id, sequences)
        _, _, lock_path = self._cache_paths(cache_key)
        with self._exclusive_file_lock(lock_path):
            cache_object = self._load_valid_cache(cache_key, multi_id)
            cache_hit = cache_object is not None
            if cache_object is None:
                self._discard_cache_entry(cache_key)
                consensus_record = self._compute_consensus_record(
                    multi_id, sequences
                )
                cache_object = self._publish_cache_record(
                    cache_key, multi_id, consensus_record
                )
            self._atomic_materialize_cache(cache_object, output_path)
        return output_path, cache_hit

    @staticmethod
    def _progress_interval(total_groups):
        """Choose an interval that emits at most about twenty summaries."""
        return max(1, total_groups // 20)

    def run(self):
        """Serialize one output generation and process its groups."""
        if not self.acquire_run_lock:
            return self._run_locked()
        with self._workspace_directory_lock():
            return self._run_locked()

    def _run_locked(self):
        """Build a complete staged generation before publishing outputs."""
        if not self.acquire_run_lock:
            self._validate_parent_workspace_lock()
        published_output_dir = self.output_dir
        staging_dir = Path(tempfile.mkdtemp(
            dir=published_output_dir,
            prefix=".alignment-generation-",
        ))
        self.output_dir = staging_dir
        try:
            result = self._process_groups()
            self.output_dir = published_output_dir
            if not self.acquire_run_lock:
                self._validate_parent_workspace_lock()
            self._publish_output_generation(
                staging_dir, published_output_dir
            )
            return result
        finally:
            self.output_dir = published_output_dir
            shutil.rmtree(staging_dir, ignore_errors=True)

    def _publish_output_generation(self, staging_dir, output_dir):
        """Replace published consensus files after a complete successful run."""
        staging_dir = Path(staging_dir)
        output_dir = Path(output_dir)
        backup_dir = Path(tempfile.mkdtemp(
            dir=output_dir,
            prefix=".previous-alignment-generation-",
        ))
        previous_outputs = list(output_dir.glob(
            "*_aligned_consensus.fa"
        ))
        staged_outputs = sorted(staging_dir.glob(
            "*_aligned_consensus.fa"
        ))
        remove_backup = False
        try:
            for previous_output in previous_outputs:
                backup_path = backup_dir / previous_output.name
                os.replace(previous_output, backup_path)
            for staged_output in staged_outputs:
                published_path = output_dir / staged_output.name
                os.replace(staged_output, published_path)
            self._fsync_directory(output_dir)
            remove_backup = True
        except BaseException as publication_error:
            rollback_errors = []
            for staged_output in reversed(staged_outputs):
                published_path = output_dir / staged_output.name
                if staged_output.exists() or not published_path.exists():
                    continue
                try:
                    os.replace(published_path, staged_output)
                except BaseException as rollback_error:
                    rollback_errors.append(rollback_error)
            for backup_path in backup_dir.glob(
                "*_aligned_consensus.fa"
            ):
                previous_output = output_dir / backup_path.name
                try:
                    os.replace(backup_path, previous_output)
                except BaseException as rollback_error:
                    rollback_errors.append(rollback_error)
            if rollback_errors:
                raise WorkspaceLockSafetyError(
                    "Alignment output publication failed and the previous "
                    "generation could not be fully restored. Preserved "
                    f"backup: {backup_dir}"
                ) from publication_error
            remove_backup = True
            raise
        finally:
            if remove_backup:
                shutil.rmtree(backup_dir, ignore_errors=True)

    def _validate_output_generation(self, expected_multi_ids):
        """Require an exact set of non-empty single-record consensus FASTAs."""
        suffix = "_aligned_consensus.fa"
        output_paths = sorted(self.output_dir.glob(f"*{suffix}"))
        observed_multi_ids = {
            output_path.name[:-len(suffix)]
            for output_path in output_paths
        }
        expected_multi_ids = set(expected_multi_ids)
        if observed_multi_ids != expected_multi_ids:
            missing = sorted(expected_multi_ids - observed_multi_ids)
            unexpected = sorted(observed_multi_ids - expected_multi_ids)
            details = []
            if missing:
                details.append("missing=" + ",".join(missing))
            if unexpected:
                details.append("unexpected=" + ",".join(unexpected))
            raise RuntimeError(
                "Alignment output generation is incomplete: "
                + "; ".join(details)
            )

        for output_path in output_paths:
            expected_multi_id = output_path.name[:-len(suffix)]
            with output_path.open("r", encoding="utf-8") as handle:
                records = list(SeqIO.parse(handle, "fasta"))
            if (
                len(records) != 1
                or records[0].id != expected_multi_id
                or len(records[0].seq) == 0
            ):
                raise RuntimeError(
                    f"Invalid aligned consensus FASTA: {output_path}"
                )

    def _process_groups(self):
        """Process all groups, retaining successful cache entries on failure."""
        print("Discovering Multi-gene groups...", flush=True)
        multi_groups = self.group_sequences_by_multi_id()
        if not multi_groups:
            print("No Multi-gene groups were found.", flush=True)
            return None

        total_groups = len(multi_groups)
        max_workers = min(self.workers, total_groups)
        print(
            f"Processing {total_groups} Multi-gene groups with "
            f"{max_workers} worker(s)...",
            flush=True,
        )

        successful_count = 0
        cache_hit_count = 0
        failures = []
        completed_count = 0
        progress_interval = self._progress_interval(total_groups)
        group_iterator = iter(multi_groups.items())

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            pending = {}

            def submit_next_group():
                try:
                    multi_id, group_source = next(group_iterator)
                except StopIteration:
                    return False
                future = executor.submit(
                    self._process_group, multi_id, group_source
                )
                pending[future] = multi_id
                return True

            for _ in range(max_workers):
                if not submit_next_group():
                    break

            while pending:
                completed_futures, _ = wait(
                    tuple(pending), return_when=FIRST_COMPLETED
                )
                for future in completed_futures:
                    multi_id = pending.pop(future)
                    try:
                        _, cache_hit = future.result()
                        successful_count += 1
                        cache_hit_count += int(cache_hit)
                    except Exception as error:
                        failures.append((multi_id, str(error)))

                    completed_count += 1
                    if (
                        completed_count % progress_interval == 0
                        or completed_count == total_groups
                    ):
                        print(
                            "Alignment progress: "
                            f"{completed_count}/{total_groups} complete, "
                            f"{successful_count} successful, "
                            f"{len(failures)} failed, "
                            f"{cache_hit_count} cache hits",
                            flush=True,
                        )
                    submit_next_group()

        print(
            "Alignment summary: "
            f"{total_groups} groups, {successful_count} successful, "
            f"{len(failures)} failed, {cache_hit_count} cache hits",
            flush=True,
        )

        if failures:
            failure_details = "; ".join(
                f"{multi_id}: {message}"
                for multi_id, message in sorted(failures)
            )
            raise RuntimeError(
                f"Alignment failed for {len(failures)} of "
                f"{total_groups} groups: {failure_details}"
            )

        self._validate_output_generation(multi_groups)

        return {
            "total_groups": total_groups,
            "successful_groups": successful_count,
            "failed_groups": 0,
        }


def _positive_integer(value):
    """Parse a positive integer for argparse."""
    try:
        parsed_value = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "workers must be a positive integer"
        ) from error
    if parsed_value < 1:
        raise argparse.ArgumentTypeError(
            "workers must be a positive integer"
        )
    return parsed_value


def main():
    """Run the alignment command-line interface."""
    parser = argparse.ArgumentParser(
        description=(
            "Align Multi-gene transcripts and generate consensus sequences."
        )
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input directory path."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output directory path."
    )
    parser.add_argument(
        "--cache-dir",
        help=(
            "Persistent cache directory. Defaults to a sibling of the output "
            "directory."
        ),
    )
    parser.add_argument(
        "--workers",
        type=_positive_integer,
        default=1,
        help="Maximum concurrent groups (default: 1).",
    )
    parser.add_argument(
        "--tool-identity",
        default=None,
        help=(
            "Tool or container identity included in cache keys. Defaults to "
            "an identity derived from this script, Python, Biopython, and "
            "MAFFT."
        ),
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Disable group cache reads and writes.",
    )
    parser.add_argument(
        "--run-lock-path",
        help=(
            "Workspace generation lock path. Defaults to "
            "<output-parent>/.annotation-workspace.lock.d."
        ),
    )
    parser.add_argument(
        "--no-run-lock",
        action="store_true",
        help="Skip locking only when a parent process holds the workspace lock.",
    )
    parser.add_argument(
        "--parent-lock-token",
        help="Ownership token supplied by the parent holding --no-run-lock.",
    )
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(
            f"Input directory does not exist: {args.input}",
            file=sys.stderr,
        )
        sys.exit(1)

    processor = AlignmentProcessor(
        args.input,
        args.output,
        cache_dir=args.cache_dir,
        workers=args.workers,
        tool_identity=args.tool_identity,
        cache_enabled=not args.no_cache,
        run_lock_path=args.run_lock_path,
        acquire_run_lock=not args.no_run_lock,
        parent_lock_token=args.parent_lock_token,
    )
    try:
        processor.run()
    except Exception as error:
        print(f"Alignment processing failed: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
