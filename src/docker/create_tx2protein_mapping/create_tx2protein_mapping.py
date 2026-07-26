#!/usr/bin/env python3
"""Create a validated tWhEEL tx2protein.tsv from RefSeq or Ensembl GFF3."""

import argparse
import gzip
import json
import os
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


AUDIT_SCHEMA_VERSION = 1
VERSION_SUFFIX_RE = re.compile(r"\.\d+$")
REFSEQ_CDS_PROTEIN_ID_RE = re.compile(r"\[protein_id=([^\]]*)\]")
TRANSCRIPT_PREFIXES = ("transcript:", "rna-", "mRNA:", "RNA:")
PROTEIN_PREFIXES = ("protein:", "cds-", "CDS:")


class MappingBuildError(ValueError):
    """Raised when the requested mapping cannot be produced safely."""


def _open_text(path):
    path = os.fspath(path)
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _strip_known_prefixes(value, prefixes):
    value = unquote(value).strip().strip('"')
    changed = True
    while changed:
        changed = False
        for prefix in prefixes:
            if value.startswith(prefix):
                value = value[len(prefix):]
                changed = True
                break
    return value.strip()


def _clean_transcript_id(value):
    return _strip_known_prefixes(value, TRANSCRIPT_PREFIXES)


def _clean_protein_id(value):
    return _strip_known_prefixes(value, PROTEIN_PREFIXES)


def _parse_attributes(text):
    attributes = defaultdict(list)
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        key, separator, raw_value = item.partition("=")
        if not separator:
            continue
        key = key.strip()
        for value in raw_value.split(","):
            value = unquote(value).strip().strip('"')
            if value:
                attributes[key].append(value)
    return dict(attributes)


def _read_fasta_headers(path, source_name):
    """Return ordered FASTA IDs and full headers after basic validation."""
    records = {}
    current_id = None
    current_line = None
    current_has_sequence = False

    def store_current():
        if current_id is None:
            return
        if not current_has_sequence:
            raise MappingBuildError(
                f"{source_name.capitalize()} FASTA record {current_id!r} has "
                f"no sequence data ({path}:{current_line})"
            )

    try:
        with _open_text(path) as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    store_current()
                    header = line[1:].strip()
                    if not header:
                        raise MappingBuildError(
                            f"Empty {source_name} FASTA header at "
                            f"{path}:{line_number}"
                        )
                    record_id = header.split(maxsplit=1)[0]
                    if record_id in records:
                        raise MappingBuildError(
                            f"Duplicate {source_name} FASTA identifier "
                            f"{record_id!r} at {path}:{line_number}"
                        )
                    records[record_id] = header
                    current_id = record_id
                    current_line = line_number
                    current_has_sequence = False
                else:
                    if current_id is None:
                        raise MappingBuildError(
                            f"Sequence data precedes the first {source_name} "
                            f"FASTA header at {path}:{line_number}"
                        )
                    current_has_sequence = True
        store_current()
    except OSError as error:
        raise MappingBuildError(
            f"Failed to read {source_name} FASTA {path}: {error}"
        ) from error

    if not records:
        raise MappingBuildError(
            f"{source_name.capitalize()} FASTA contains no records: {path}"
        )
    return records


class IdentifierResolver:
    """Resolve exact IDs first, then a unique version-insensitive match."""

    def __init__(self, identifiers):
        self.identifiers = set(identifiers)
        self.by_versionless = defaultdict(list)
        for identifier in identifiers:
            self.by_versionless[VERSION_SUFFIX_RE.sub("", identifier)].append(
                identifier
            )

    def resolve(self, candidate):
        if candidate in self.identifiers:
            return candidate, "exact"
        versionless = VERSION_SUFFIX_RE.sub("", candidate)
        matches = self.by_versionless.get(versionless, [])
        if len(matches) == 1:
            match = matches[0]
            candidate_has_version = versionless != candidate
            match_has_version = VERSION_SUFFIX_RE.sub("", match) != match
            if candidate_has_version and match_has_version:
                return None, "conflicting-version"
            return match, "version-normalized"
        if len(matches) > 1:
            return None, "ambiguous-version"
        return None, "absent"


def _extract_gff_candidates(attributes):
    proteins = list(attributes.get("protein_id", []))
    transcripts = list(attributes.get("transcript_id", []))
    if not transcripts:
        transcripts = list(attributes.get("Parent", []))

    cleaned_proteins = {
        _clean_protein_id(value) for value in proteins if value.strip()
    }
    cleaned_transcripts = {
        _clean_transcript_id(value) for value in transcripts if value.strip()
    }
    cleaned_proteins.discard("")
    cleaned_transcripts.discard("")
    return cleaned_transcripts, cleaned_proteins


def _detect_source(scores):
    if scores["refseq"] <= 0 and scores["ensembl"] <= 0:
        return None
    if scores["refseq"] == scores["ensembl"]:
        return None
    return max(scores, key=scores.get)


def _parse_gff(path):
    raw_pairs = set()
    scores = {"refseq": 0, "ensembl": 0}
    counts = {
        "gff_feature_rows": 0,
        "gff_cds_rows": 0,
        "cds_rows_without_protein_id_attribute": 0,
        "cds_rows_without_protein_candidate": 0,
        "cds_rows_without_transcript_candidate": 0,
    }

    try:
        with _open_text(path) as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("#"):
                    lower = line.lower()
                    if "annotation-source ncbi refseq" in lower:
                        scores["refseq"] += 100
                    elif "processor ncbi" in lower:
                        scores["refseq"] += 30
                    continue

                fields = line.split("\t")
                if len(fields) != 9:
                    raise MappingBuildError(
                        f"Malformed GFF3 row at {path}:{line_number}: "
                        f"expected 9 tab-separated columns, found {len(fields)}"
                    )
                counts["gff_feature_rows"] += 1
                if fields[2].lower() != "cds":
                    continue
                counts["gff_cds_rows"] += 1
                attributes = _parse_attributes(fields[8])

                parent_values = attributes.get("Parent", [])
                id_values = attributes.get("ID", [])
                if any(value.startswith("rna-") for value in parent_values):
                    scores["refseq"] += 3
                if any(value.startswith("transcript:") for value in parent_values):
                    scores["ensembl"] += 3
                if any(value.startswith("cds-") for value in id_values):
                    scores["refseq"] += 1
                if any(value.startswith("CDS:") for value in id_values):
                    scores["ensembl"] += 1
                if "gbkey" in attributes:
                    scores["refseq"] += 1
                if fields[1] in {"RefSeq", "BestRefSeq", "Gnomon"}:
                    scores["refseq"] += 1

                if not attributes.get("protein_id"):
                    counts["cds_rows_without_protein_id_attribute"] += 1
                transcripts, proteins = _extract_gff_candidates(attributes)
                if not proteins:
                    counts["cds_rows_without_protein_candidate"] += 1
                    continue
                if not transcripts:
                    counts["cds_rows_without_transcript_candidate"] += 1
                    continue
                for transcript_id in transcripts:
                    for protein_id in proteins:
                        raw_pairs.add((transcript_id, protein_id))
    except OSError as error:
        raise MappingBuildError(f"Failed to read GFF3 {path}: {error}") from error

    if counts["gff_cds_rows"] == 0:
        raise MappingBuildError(f"GFF3 contains no CDS rows: {path}")
    counts["unique_gff_candidate_pairs"] = len(raw_pairs)
    return raw_pairs, scores, counts


def _ensembl_header_pairs(protein_headers):
    pairs = set()
    for protein_id, header in protein_headers.items():
        for token in header.split()[1:]:
            if token.startswith("transcript:"):
                transcript_id = _clean_transcript_id(token)
                if transcript_id:
                    pairs.add((transcript_id, protein_id))
    return pairs


def _refseq_cds_header_pairs(transcript_headers, max_examples):
    pairs = set()
    missing_protein_ids = []
    multiple_protein_ids = []
    empty_protein_ids = []
    for transcript_id, header in transcript_headers.items():
        raw_protein_ids = REFSEQ_CDS_PROTEIN_ID_RE.findall(header)
        if not raw_protein_ids:
            missing_protein_ids.append(transcript_id)
            continue
        if len(raw_protein_ids) > 1:
            multiple_protein_ids.append(
                f"{transcript_id}->{'|'.join(raw_protein_ids)}"
            )
            continue
        protein_id = _clean_protein_id(raw_protein_ids[0])
        if not protein_id:
            empty_protein_ids.append(transcript_id)
            continue
        pairs.add((transcript_id, protein_id))

    counts = {
        "refseq_cds_fasta_candidate_pairs": len(pairs),
        "refseq_cds_fasta_records_without_protein_id": len(
            missing_protein_ids
        ),
        "refseq_cds_fasta_records_with_multiple_protein_ids": len(
            multiple_protein_ids
        ),
        "refseq_cds_fasta_records_with_empty_protein_id": len(
            empty_protein_ids
        ),
    }
    examples = {
        "refseq_cds_fasta_records_without_protein_id": _bounded(
            missing_protein_ids, max_examples
        ),
        "refseq_cds_fasta_records_with_multiple_protein_ids": _bounded(
            multiple_protein_ids, max_examples
        ),
        "refseq_cds_fasta_records_with_empty_protein_id": _bounded(
            empty_protein_ids, max_examples
        ),
    }
    problems = []
    if multiple_protein_ids:
        problems.append(
            _problem(
                "multiple_protein_ids_in_refseq_cds_fasta_header",
                multiple_protein_ids,
                max_examples,
            )
        )
    if empty_protein_ids:
        problems.append(
            _problem(
                "empty_protein_id_in_refseq_cds_fasta_header",
                empty_protein_ids,
                max_examples,
            )
        )
    return pairs, counts, examples, problems


def _bounded(values, limit):
    values = list(values)
    return values[:limit]


def _problem(code, values, limit):
    values = list(values)
    return {
        "code": code,
        "count": len(values),
        "examples": _bounded(values, limit),
    }


def _make_audit(args, audit_path, filtered_protein_path):
    return {
        "schema_version": AUDIT_SCHEMA_VERSION,
        "status": "running",
        "source_requested": args.source,
        "source_detected": None,
        "refseq_xp": args.refseq_xp,
        "refseq_transcript_source": args.refseq_transcript_source,
        "unresolved_transcripts": args.unresolved_transcripts,
        "unmapped_proteins": args.unmapped_proteins,
        "inputs": {
            "gff3": os.path.abspath(args.annotation),
            "transcripts": os.path.abspath(args.transcripts),
            "proteins": os.path.abspath(args.proteins),
        },
        "outputs": {
            "mapping": os.path.abspath(args.output),
            "audit": os.path.abspath(audit_path),
            "filtered_proteins": (
                os.path.abspath(filtered_protein_path)
                if filtered_protein_path
                else None
            ),
        },
        "counts": {},
        "problems": [],
    }


def _build_mapping(args, audit):
    transcript_headers = _read_fasta_headers(args.transcripts, "transcript")
    protein_headers = _read_fasta_headers(args.proteins, "protein")
    raw_pairs, scores, gff_counts = _parse_gff(args.annotation)
    detected_source = _detect_source(scores)
    audit["source_detected"] = detected_source
    audit["source_scores"] = scores
    audit["counts"].update(gff_counts)
    audit["counts"].update(
        {
            "transcript_fasta_records": len(transcript_headers),
            "protein_fasta_records": len(protein_headers),
        }
    )

    if detected_source is None:
        audit["problems"].append(
            {
                "code": "source_detection_failed",
                "count": 1,
                "examples": [str(scores)],
            }
        )
        return None, None
    if args.source != "auto" and args.source != detected_source:
        audit["problems"].append(
            {
                "code": "source_mismatch",
                "count": 1,
                "examples": [
                    f"requested={args.source}, detected={detected_source}"
                ],
            }
        )
        return None, None
    source = detected_source if args.source == "auto" else args.source
    if args.refseq_xp == "exclude" and source != "refseq":
        audit["problems"].append(
            {
                "code": "refseq_xp_option_requires_refseq",
                "count": 1,
                "examples": [f"detected source: {source}"],
            }
        )
        return None, None
    if args.refseq_transcript_source == "cds-fasta" and source != "refseq":
        audit["problems"].append(
            {
                "code": "refseq_cds_fasta_option_requires_refseq",
                "count": 1,
                "examples": [f"detected source: {source}"],
            }
        )
        return None, None

    header_examples = {}
    if args.refseq_transcript_source == "cds-fasta":
        (
            raw_pairs,
            header_counts,
            header_examples,
            header_problems,
        ) = _refseq_cds_header_pairs(
            transcript_headers, args.max_examples
        )
        audit["counts"].update(header_counts)
        audit["problems"].extend(header_problems)
    elif source == "ensembl":
        header_pairs = _ensembl_header_pairs(protein_headers)
        raw_pairs.update(header_pairs)
        audit["counts"]["ensembl_protein_header_pairs"] = len(header_pairs)

    all_protein_ids = set(protein_headers)
    excluded_xp = sorted(
        protein_id for protein_id in all_protein_ids
        if args.refseq_xp == "exclude" and protein_id.startswith("XP_")
    )
    selected_protein_ids = all_protein_ids - set(excluded_xp)
    audit["counts"]["excluded_xp_proteins"] = len(excluded_xp)
    audit["counts"]["protein_records_after_xp_policy"] = len(
        selected_protein_ids
    )
    if not selected_protein_ids:
        audit["problems"].append(
            {
                "code": "no_selected_proteins",
                "count": 1,
                "examples": ["No protein remains after applying the XP policy"],
            }
        )
        return None, None

    protein_resolver = IdentifierResolver(protein_headers)
    transcript_resolver = IdentifierResolver(transcript_headers)
    mappings = defaultdict(set)
    absent_annotation_proteins = set()
    ambiguous_protein_identifiers = set()
    ambiguous_transcript_identifiers = set()
    conflicting_protein_versions = set()
    conflicting_transcript_versions = set()
    unresolved_transcripts = set()
    pair_modes = {}

    for raw_transcript, raw_protein in raw_pairs:
        protein_id, protein_mode = protein_resolver.resolve(raw_protein)
        if protein_id is None:
            if protein_mode == "ambiguous-version":
                ambiguous_protein_identifiers.add(raw_protein)
            elif protein_mode == "conflicting-version":
                conflicting_protein_versions.add(raw_protein)
            else:
                absent_annotation_proteins.add(raw_protein)
            continue
        if protein_id not in selected_protein_ids:
            continue
        transcript_id, transcript_mode = transcript_resolver.resolve(
            raw_transcript
        )
        if transcript_id is None:
            unresolved = (protein_id, raw_transcript, transcript_mode)
            if transcript_mode == "ambiguous-version":
                ambiguous_transcript_identifiers.add(unresolved)
            elif transcript_mode == "conflicting-version":
                conflicting_transcript_versions.add(unresolved)
            else:
                unresolved_transcripts.add(unresolved)
            continue
        mappings[protein_id].add(transcript_id)
        pair = (transcript_id, protein_id)
        pair_mode = (
            "exact"
            if protein_mode == "exact" and transcript_mode == "exact"
            else "version-normalized"
        )
        if pair_mode == "exact" or pair not in pair_modes:
            pair_modes[pair] = pair_mode

    missing_proteins = sorted(selected_protein_ids - set(mappings))
    ambiguous = sorted(
        f"{protein_id}->{'|'.join(sorted(transcript_ids))}"
        for protein_id, transcript_ids in mappings.items()
        if len(transcript_ids) != 1
    )
    excluded_unmapped = []
    if missing_proteins and args.unmapped_proteins == "exclude":
        excluded_unmapped = missing_proteins
        selected_protein_ids -= set(excluded_unmapped)
        missing_proteins = []
    if not selected_protein_ids:
        audit["problems"].append(
            {
                "code": "no_selected_proteins",
                "count": 1,
                "examples": [
                    "No protein remains after applying the XP and unmapped "
                    "protein policies"
                ],
            }
        )
    if missing_proteins:
        audit["problems"].append(
            _problem("unmapped_proteins", missing_proteins, args.max_examples)
        )
    if ambiguous:
        audit["problems"].append(
            _problem("ambiguous_protein_mappings", ambiguous, args.max_examples)
        )
    if ambiguous_protein_identifiers:
        audit["problems"].append(
            _problem(
                "ambiguous_protein_identifiers",
                sorted(ambiguous_protein_identifiers),
                args.max_examples,
            )
        )
    if ambiguous_transcript_identifiers:
        values = sorted(
            f"{protein_id}->{raw_transcript} ({mode})"
            for protein_id, raw_transcript, mode
            in ambiguous_transcript_identifiers
        )
        audit["problems"].append(
            _problem(
                "ambiguous_transcript_identifiers",
                values,
                args.max_examples,
            )
        )
    if conflicting_protein_versions:
        audit["problems"].append(
            _problem(
                "conflicting_protein_versions",
                sorted(conflicting_protein_versions),
                args.max_examples,
            )
        )
    if conflicting_transcript_versions:
        values = sorted(
            f"{protein_id}->{raw_transcript} ({mode})"
            for protein_id, raw_transcript, mode
            in conflicting_transcript_versions
        )
        audit["problems"].append(
            _problem(
                "conflicting_transcript_versions",
                values,
                args.max_examples,
            )
        )
    unresolved_selected = sorted(
        f"{protein_id}->{raw_transcript} ({mode})"
        for protein_id, raw_transcript, mode in unresolved_transcripts
        if protein_id in selected_protein_ids
    )
    excluded_unresolved_transcripts = []
    if args.unresolved_transcripts == "exclude":
        excluded_unresolved_transcripts = sorted(
            f"{protein_id}->{raw_transcript} ({mode})"
            for protein_id, raw_transcript, mode in unresolved_transcripts
        )
    elif unresolved_selected:
        audit["problems"].append(
            _problem(
                "unresolved_transcript_ids",
                unresolved_selected,
                args.max_examples,
            )
        )

    audit["counts"].update(
        {
            "annotation_protein_ids_absent_from_input": len(
                absent_annotation_proteins
            ),
            "exact_resolved_pairs": sum(
                mode == "exact" for mode in pair_modes.values()
            ),
            "version_normalized_pairs": sum(
                mode == "version-normalized" for mode in pair_modes.values()
            ),
            "mapped_proteins": len(mappings),
            "excluded_unmapped_proteins": len(excluded_unmapped),
            "excluded_unresolved_transcript_pairs": len(
                excluded_unresolved_transcripts
            ),
            "selected_protein_records": len(selected_protein_ids),
            "mapped_transcripts": len(
                {transcript for values in mappings.values() for transcript in values}
            ),
        }
    )
    audit["examples"] = {
        "excluded_xp_proteins": _bounded(excluded_xp, args.max_examples),
        "excluded_unmapped_proteins": _bounded(
            excluded_unmapped, args.max_examples
        ),
        "excluded_unresolved_transcript_pairs": _bounded(
            excluded_unresolved_transcripts, args.max_examples
        ),
        "annotation_protein_ids_absent_from_input": _bounded(
            sorted(absent_annotation_proteins), args.max_examples
        ),
    }
    audit["examples"].update(header_examples)
    if audit["problems"]:
        return None, selected_protein_ids

    rows = sorted(
        (next(iter(transcript_ids)), protein_id)
        for protein_id, transcript_ids in mappings.items()
    )
    if not rows:
        audit["problems"].append(
            {
                "code": "no_output_mappings",
                "count": 1,
                "examples": ["No transcript-protein mapping can be emitted"],
            }
        )
        return None, selected_protein_ids
    return rows, selected_protein_ids


def _temporary_path(target):
    target = Path(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=str(target.parent),
        prefix=f".{target.name}.",
        suffix=".tmp",
        delete=False,
    )
    return handle


def _write_mapping_temp(rows, target):
    handle = _temporary_path(target)
    try:
        with handle:
            for transcript_id, protein_id in rows:
                handle.write(f"{transcript_id}\t{protein_id}\n")
        return handle.name
    except Exception:
        try:
            os.unlink(handle.name)
        except OSError:
            pass
        raise


def _write_filtered_proteins_temp(source, selected_ids, target):
    handle = _temporary_path(target)
    written_ids = set()
    include_record = False
    try:
        with handle, _open_text(source) as input_handle:
            for raw_line in input_handle:
                if raw_line.startswith(">"):
                    header = raw_line[1:].strip()
                    record_id = header.split(maxsplit=1)[0]
                    include_record = record_id in selected_ids
                    if include_record:
                        written_ids.add(record_id)
                if include_record:
                    output_line = (
                        raw_line if raw_line.endswith("\n") else raw_line + "\n"
                    )
                    handle.write(output_line)
        missing = sorted(selected_ids - written_ids)
        if missing:
            raise MappingBuildError(
                "Failed to copy selected protein FASTA records: "
                + ", ".join(missing[:5])
            )
        return handle.name
    except Exception:
        try:
            os.unlink(handle.name)
        except OSError:
            pass
        raise


def _write_json_atomic(data, target):
    handle = _temporary_path(target)
    try:
        with handle:
            json.dump(data, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(handle.name, target)
    except Exception:
        try:
            os.unlink(handle.name)
        except OSError:
            pass
        raise


def _resolved(path):
    return str(Path(path).expanduser().resolve())


def _check_paths(args, audit_path, filtered_protein_path):
    inputs = {
        _resolved(args.annotation),
        _resolved(args.transcripts),
        _resolved(args.proteins),
    }
    outputs = [_resolved(args.output), _resolved(audit_path)]
    if filtered_protein_path:
        outputs.append(_resolved(filtered_protein_path))
    if len(set(outputs)) != len(outputs):
        raise MappingBuildError("Output and audit paths must be distinct")
    collisions = sorted(inputs & set(outputs))
    if collisions:
        raise MappingBuildError(
            "Output paths must not overwrite input files: " + ", ".join(collisions)
        )
    if not args.force:
        existing = [path for path in outputs if os.path.exists(path)]
        if existing:
            raise MappingBuildError(
                "Refusing to overwrite existing output(s) without --force: "
                + ", ".join(existing)
            )


def _default_filtered_protein_path(output_path):
    output = Path(output_path)
    stem = output.stem if output.suffix else output.name
    return str(output.with_name(f"{stem}.matched.faa"))


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Create a headerless transcript_id<TAB>protein_id mapping from "
            "matching RefSeq or Ensembl annotation and sequence files."
        )
    )
    parser.add_argument(
        "--annotation",
        "--gff",
        required=True,
        help="GFF3 annotation file (.gz is supported)",
    )
    parser.add_argument(
        "--transcripts",
        required=True,
        help=(
            "Transcript/cDNA FASTA, or an NCBI cds_from_genomic nucleotide "
            "FASTA in RefSeq CDS-header mode (.gz is supported)"
        ),
    )
    parser.add_argument(
        "--proteins",
        required=True,
        help="Protein/peptide FASTA (.gz is supported)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Headerless tx2protein.tsv output",
    )
    parser.add_argument(
        "--source",
        choices=("auto", "refseq", "ensembl"),
        default="auto",
        help="Annotation source (default: detect RefSeq or Ensembl)",
    )
    parser.add_argument(
        "--refseq-xp",
        choices=("include", "exclude"),
        default="include",
        help=(
            "Include XP_ proteins, or exclude them and create a matching "
            "filtered protein FASTA (default: include)"
        ),
    )
    parser.add_argument(
        "--refseq-transcript-source",
        choices=("gff3", "cds-fasta"),
        default="gff3",
        help=(
            "Obtain RefSeq transcript-protein links from GFF3 CDS attributes, "
            "or from [protein_id=...] fields in an NCBI cds_from_genomic "
            "FASTA supplied with --transcripts (default: gff3)"
        ),
    )
    parser.add_argument(
        "--unresolved-transcripts",
        choices=("error", "exclude"),
        default="error",
        help=(
            "Fail on GFF3 transcript candidates absent from the transcript "
            "FASTA, or explicitly exclude only those candidates while "
            "retaining uniquely resolved mappings (default: error)"
        ),
    )
    parser.add_argument(
        "--filtered-proteins",
        help=(
            "Output path for the protein FASTA matching the generated map; "
            "valid with --refseq-xp exclude or --unmapped-proteins exclude "
            "(default: <mapping stem>.matched.faa)"
        ),
    )
    parser.add_argument(
        "--unmapped-proteins",
        choices=("error", "exclude"),
        default="error",
        help=(
            "Fail on proteins without a resolvable transcript, or exclude "
            "them and create a matching filtered FASTA (default: error)"
        ),
    )
    parser.add_argument(
        "--audit",
        help="Audit JSON output (default: <output>.audit.json)",
    )
    parser.add_argument(
        "--max-examples",
        type=int,
        default=20,
        help="Maximum identifiers retained per audit example list (default: 20)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Atomically replace existing output files after successful validation",
    )
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.max_examples < 1:
        parser.error("--max-examples must be at least 1")
    filtering_requested = (
        args.refseq_xp == "exclude" or args.unmapped_proteins == "exclude"
    )
    if args.filtered_proteins and not filtering_requested:
        parser.error(
            "--filtered-proteins requires --refseq-xp exclude or "
            "--unmapped-proteins exclude"
        )

    audit_path = args.audit or f"{args.output}.audit.json"
    filtered_protein_path = None
    if filtering_requested:
        filtered_protein_path = (
            args.filtered_proteins or _default_filtered_protein_path(args.output)
        )

    try:
        _check_paths(args, audit_path, filtered_protein_path)
    except MappingBuildError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2

    audit = _make_audit(args, audit_path, filtered_protein_path)
    mapping_temp = None
    protein_temp = None
    try:
        rows, selected_protein_ids = _build_mapping(args, audit)
        if audit["problems"]:
            audit["status"] = "error"
            _write_json_atomic(audit, audit_path)
            print(
                f"ERROR: mapping validation failed; no mapping was written. "
                f"See {audit_path}",
                file=sys.stderr,
            )
            return 2

        mapping_temp = _write_mapping_temp(rows, args.output)
        if filtered_protein_path:
            protein_temp = _write_filtered_proteins_temp(
                args.proteins, selected_protein_ids, filtered_protein_path
            )
        os.replace(mapping_temp, args.output)
        mapping_temp = None
        if protein_temp:
            os.replace(protein_temp, filtered_protein_path)
            protein_temp = None
        audit["status"] = "ok"
        audit["counts"]["output_mapping_rows"] = len(rows)
        _write_json_atomic(audit, audit_path)
    except (MappingBuildError, OSError, UnicodeError) as error:
        audit["status"] = "error"
        audit["problems"].append(
            {"code": "runtime_error", "count": 1, "examples": [str(error)]}
        )
        try:
            _write_json_atomic(audit, audit_path)
        except OSError:
            pass
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    finally:
        for temporary in (mapping_temp, protein_temp):
            if temporary:
                try:
                    os.unlink(temporary)
                except OSError:
                    pass

    print(f"Source: {audit['source_detected']}")
    print(f"Mappings: {audit['counts']['output_mapping_rows']}")
    print(f"Mapping: {args.output}")
    if filtered_protein_path:
        print(f"Matched proteins: {filtered_protein_path}")
    print(f"Audit: {audit_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
