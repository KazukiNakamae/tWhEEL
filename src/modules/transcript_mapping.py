"""Protein-to-transcript mapping loading and validation."""

import os


class TranscriptMappingError(ValueError):
    """Raised when protein and transcript identifiers are inconsistent."""


def _format_examples(values, limit=5):
    """Format a bounded, deterministic identifier sample for an error."""
    values = list(values)
    displayed = values[:limit]
    suffix = '' if len(values) <= limit else f", ... (+{len(values) - limit})"
    return ', '.join(displayed) + suffix


def load_mapping_file(mapping_path):
    """Load a headerless transcript_id-to-protein_id TSV as protein-to-transcript."""
    mapping_path = os.fspath(mapping_path)
    mapping = {}
    source_lines = {}

    try:
        with open(mapping_path, 'r', encoding='utf-8') as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.strip()
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 2:
                    raise TranscriptMappingError(
                        f"Invalid protein-transcript mapping at "
                        f"{mapping_path}:{line_number}: expected a headerless "
                        "transcript_id<TAB>protein_id row"
                    )

                transcript_id = fields[0].strip()
                protein_id = fields[1].strip()
                if not transcript_id or not protein_id:
                    raise TranscriptMappingError(
                        f"Invalid protein-transcript mapping at "
                        f"{mapping_path}:{line_number}: transcript_id and "
                        "protein_id must both be non-empty"
                    )

                if (
                    transcript_id.lower() == 'transcript_id'
                    and protein_id.lower() == 'protein_id'
                ):
                    raise TranscriptMappingError(
                        f"Invalid protein-transcript mapping at "
                        f"{mapping_path}:{line_number}: header rows are not "
                        "allowed; expected transcript_id<TAB>protein_id data"
                    )

                previous_transcript = mapping.get(protein_id)
                if (
                    previous_transcript is not None
                    and previous_transcript != transcript_id
                ):
                    raise TranscriptMappingError(
                        "Conflicting protein-transcript mappings for "
                        f"protein_id {protein_id!r}: {previous_transcript!r} "
                        f"at line {source_lines[protein_id]} and "
                        f"{transcript_id!r} at line {line_number}"
                    )

                mapping[protein_id] = transcript_id
                source_lines.setdefault(protein_id, line_number)
    except OSError as error:
        raise RuntimeError(
            f"Failed to read protein-transcript mapping: {mapping_path}: "
            f"{error}"
        ) from error

    if not mapping:
        raise TranscriptMappingError(
            "Protein-transcript mapping contains no data rows: "
            f"{mapping_path}"
        )

    return mapping


def _read_fasta_data(fasta_path, source_name, include_sequences):
    """Read validated FASTA data without retaining sequences when unnecessary."""
    fasta_path = os.fspath(fasta_path)
    records = {}
    current_id = None
    current_value = [] if include_sequences else 0
    current_line = None

    def store_current_record():
        if current_id is None:
            return
        if current_id in records:
            raise TranscriptMappingError(
                f"Duplicate {source_name} FASTA identifier {current_id!r} "
                f"at {fasta_path}:{current_line}"
            )
        records[current_id] = (
            ''.join(current_value) if include_sequences else current_value
        )

    try:
        with open(fasta_path, 'r', encoding='utf-8') as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    store_current_record()
                    header = line[1:].strip()
                    if not header:
                        raise TranscriptMappingError(
                            f"Empty {source_name} FASTA header at "
                            f"{fasta_path}:{line_number}"
                        )
                    current_id = header.split(maxsplit=1)[0]
                    current_value = [] if include_sequences else 0
                    current_line = line_number
                else:
                    if current_id is None:
                        raise TranscriptMappingError(
                            f"Sequence data precedes the first {source_name} "
                            f"FASTA header at {fasta_path}:{line_number}"
                        )
                    if include_sequences:
                        current_value.append(line)
                    else:
                        current_value += len(line)
        store_current_record()
    except OSError as error:
        raise RuntimeError(
            f"Failed to read {source_name} FASTA: {fasta_path}: {error}"
        ) from error

    if not records:
        raise TranscriptMappingError(
            f"{source_name.capitalize()} FASTA contains no records: "
            f"{fasta_path}"
        )

    empty_ids = sorted(
        record_id for record_id, value in records.items() if not value
    )
    if empty_ids:
        raise TranscriptMappingError(
            f"{source_name.capitalize()} FASTA contains empty sequence "
            f"record(s): {_format_examples(empty_ids)}"
        )

    return records


def read_fasta_sequences(fasta_path, source_name):
    """Return validated FASTA sequences keyed by the first header token."""
    return _read_fasta_data(fasta_path, source_name, include_sequences=True)


def read_fasta_records(fasta_path, source_name):
    """Return validated FASTA record lengths keyed by identifier."""
    return _read_fasta_data(fasta_path, source_name, include_sequences=False)


def validate_mapping_against_fastas(
    mapping,
    protein_fasta_path,
    transcript_fasta_path,
):
    """Require exact protein coverage and resolvable transcript targets."""
    protein_records = read_fasta_records(protein_fasta_path, 'protein')
    transcript_records = read_fasta_records(
        transcript_fasta_path, 'transcript'
    )
    protein_ids = set(protein_records)
    mapping_protein_ids = set(mapping)

    missing_proteins = sorted(protein_ids - mapping_protein_ids)
    unexpected_proteins = sorted(mapping_protein_ids - protein_ids)
    missing_transcripts = sorted(
        (protein_id, mapping[protein_id])
        for protein_id in protein_ids & mapping_protein_ids
        if mapping[protein_id] not in transcript_records
    )

    problems = []
    if missing_proteins:
        problems.append(
            f"protein ID(s) without a mapping ({len(missing_proteins)}): "
            f"{_format_examples(missing_proteins)}"
        )
    if unexpected_proteins:
        problems.append(
            "mapping protein ID(s) absent from query.pep "
            f"({len(unexpected_proteins)}): "
            f"{_format_examples(unexpected_proteins)}"
        )
    if missing_transcripts:
        missing_pairs = [
            f"{protein_id}->{transcript_id}"
            for protein_id, transcript_id in missing_transcripts
        ]
        problems.append(
            "mapped transcript ID(s) absent from transcriptome.fa "
            f"({len(missing_pairs)}): {_format_examples(missing_pairs)}"
        )

    if problems:
        raise TranscriptMappingError(
            "Protein-transcript mapping is inconsistent with the prepared "
            "FASTA inputs: " + '; '.join(problems) + ". Expected a headerless "
            "transcript_id<TAB>protein_id mapping with exact query protein "
            "coverage. Correct the mapping and rerun the annotation step."
        )

    return {
        'protein_count': len(protein_records),
        'mapped_transcript_count': len(set(mapping.values())),
        'transcript_count': len(transcript_records),
    }
