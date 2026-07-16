#!/usr/bin/env python3
"""
Annotation module
Functional annotation via FASTA search against UniProt database
"""

import os
import subprocess
import logging
import re
from pathlib import Path
import shutil

from .fasta_search_cache import FastaSearchCache
from .input_mode import EXISTING_PROTEIN_MODE, determine_protein_input_mode
from .multi_gene_alignment import (
    MultiGeneAlignmentRunner,
    annotation_workspace_lock,
    clear_multi_gene_artifacts,
    clear_single_gene_artifacts,
)
from .transcript_mapping import (
    TranscriptMappingError,
    load_mapping_file,
    read_fasta_sequences,
    validate_mapping_against_fastas,
)


class Annotation:
    """Annotation processing class"""
    
    def __init__(self, config):
        """
        Initialize
        
        Args:
            config (dict): Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        self.protein_transcript_mapping = None
        
    def prepare_protein_db_for_goanalysis(self):
        """Prepare protein database for GO analysis"""
        protein_db_path = self.config['input']['protein_db_for_goanalysis']
        
        # Output file path
        output_path = os.path.join(self.intermediate_dir, 'protein_db.fasta')
        
        # Extract and copy the files
        if protein_db_path.endswith('.gz'):
            self.logger.info("Decompressing protein database file...")
            import gzip
            with gzip.open(protein_db_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            self.logger.info("Copying protein database file...")
            shutil.copy2(protein_db_path, output_path)
        
        self.logger.info(f"Protein database preparation complete: {output_path}")
        return output_path
    
    def load_protein_transcript_mapping(self):
        """Load protein to transcript ID mapping file"""
        mapping_path = self.config['input'].get('protein_transcript_mapping', '')
        
        if not mapping_path or not os.path.exists(mapping_path):
            self.logger.info("Protein to transcript mapping file not specified")
            return None
        
        self.logger.info(f"Loading protein to transcript mapping file: {mapping_path}")

        mapping_dict = load_mapping_file(mapping_path)
        self.logger.info(
            f"Built the mapping dictionary: {len(mapping_dict)} entries"
        )
        self.protein_transcript_mapping = mapping_dict
        return mapping_dict

    def validate_protein_transcript_mapping(
        self,
        protein_fasta_path,
        transcript_fasta_path,
    ):
        """Validate mapping coverage against the prepared pipeline FASTAs."""
        if self.protein_transcript_mapping is None:
            raise TranscriptMappingError(
                "Protein-transcript mapping has not been loaded"
            )

        validation = validate_mapping_against_fastas(
            self.protein_transcript_mapping,
            protein_fasta_path,
            transcript_fasta_path,
        )
        self.logger.info(
            "Protein-transcript mapping validation completed: "
            f"{validation['protein_count']} proteins, "
            f"{validation['mapped_transcript_count']} mapped transcripts"
        )
        return validation
    
    def run_fasta_search(self, query_path, db_path):
        """Execute FASTA search"""
        self.logger.info("Start the FASTA search against the UniProt database...")
        
        # Run ggsearch36 using intermediate_dir as the working directory
        query_filename = os.path.basename(query_path)
        db_filename = os.path.basename(db_path)
        
        cmd = [
            "docker", "run",
            "-m", "20G",
            "-v", f"{os.path.abspath(self.intermediate_dir)}:/fasta_work",
            "-w", "/fasta_work",
            "--rm",
            self.config['fasta_search']['docker_image'],
            "ggsearch36",
            "-Q", "-T", str(self.config['fasta_search']['threads']),
            "-d", "1", "-m10",
            "-E", str(self.config['fasta_search']['evalue']),
            query_filename, db_filename
        ]
        
        self.logger.debug(f"FASTA search command: {' '.join(cmd)}")

        cache = FastaSearchCache(
            self.intermediate_dir,
            namespace='annotation',
            logger=self.logger,
        )
        return cache.run(
            cmd,
            query_path,
            db_path,
            search_settings={
                'docker_image': self.config['fasta_search']['docker_image'],
                'threads': self.config['fasta_search']['threads'],
                'evalue': self.config['fasta_search']['evalue'],
                'max_hits': 1,
                'output_format': 'm10',
            },
        )
    
    def parse_fasta_results(self, search_results_path):
        """Parse FASTA search results"""
        self.logger.info("Parsing the FASTA search results...")
        
        parsed_path = os.path.join(self.intermediate_dir, 'annotation_table.txt')
        
        # 結果パース用のPerlスクリプト相当の処理をPythonで実装
        self._parse_ggsearch_output(search_results_path, parsed_path)
        
        self.logger.info(f"Save the annotation results: {parsed_path}")
        return parsed_path
    
    def _parse_ggsearch_output(self, input_path, output_path):
        """Parse ggsearch output and create annotation table.

        ggsearch36 -m10 output has two sections per query:
        1. Standard summary section with "The best scores are:" - descriptions are TRUNCATED
           (e.g., "Zinc finger protein 708 OS=H")
        2. Machine-readable (-m10) alignment section with ">>subject_id" headers - descriptions
           are FULL (e.g., "Zinc finger protein 708 OS=Homo sapiens OX=9606 GN=ZNF708")

        This method performs a two-pass parse:
        - Pass 1: Reads the standard "best scores" section to obtain subject_id, length, bits, evalue
          for the best (first) hit per query.
        - Pass 2: Reads the -m10 alignment header lines to obtain the FULL description, including
          the GN= field when available, then updates the annotations collected in Pass 1.

        Subject identifiers are treated as opaque FASTA identifiers. This supports Swiss-Prot
        (sp|), TrEMBL (tr|), NCBI pipe-delimited identifiers, and custom identifiers.
        """
        annotations = {}  # query_id -> annotation dict

        with open(input_path, 'r') as f:
            content = f.read()

        self.logger.debug(f"ggsearch output first 500 chars: {content[:500]}")

        lines = content.split('\n')

        # ------------------------------------------------------------------
        # Pass 1: Parse the standard-format section (best scores per query)
        # The query header looks like:
        #   "  1>>>NM_001001130.3.p1 GENE... - 518 aa"
        # A best-score hit line looks like:
        #   "subject_id Description ( 563) 1854 97.9 0"
        # ------------------------------------------------------------------
        current_query = None
        processing_hits = False

        for line in lines:
            line = line.strip()

            # The FASTA description is optional, so capture only the first
            # whitespace-delimited identifier before the final length field.
            query_match = re.match(
                r'^\d+>>>(?P<query_id>\S+)(?:\s+.*?)?\s+-\s+\d+\s+aa(?:\s|$)',
                line,
            )
            if query_match:
                current_query = query_match.group('query_id')
                processing_hits = False
                self.logger.debug(f"Found query: {current_query}")
                continue

            # "The best scores are:" marks the start of the hit summary for this query
            if "The best scores are:" in line:
                processing_hits = True
                continue

            # Parse the first hit without assuming a database-specific identifier prefix.
            if processing_hits and current_query and current_query not in annotations:
                hit_match = re.match(
                    r'^(?P<subject_id>\S+)(?:\s+(?P<description>.*?))?\s+'
                    r'\(\s*(?P<length>\d+)\)\s+'
                    r'(?P<raw_score>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+'
                    r'(?P<bits>[+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+'
                    r'(?P<evalue>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)$',
                    line
                )
                if hit_match:
                    subject_id = hit_match.group('subject_id')
                    description = (hit_match.group('description') or '').strip()
                    length = hit_match.group('length')
                    bits = hit_match.group('bits')
                    evalue = hit_match.group('evalue')

                    annotations[current_query] = {
                        'query_id':   current_query,
                        'subject_id': subject_id,
                        'identity':   '',
                        'length':     length,
                        'evalue':     evalue,
                        'bit_score':  bits,
                        'description': description,  # Truncated (fallback if Pass 2 fails)
                    }
                    processing_hits = False
                    self.logger.debug(f"Pass 1 - best hit: {subject_id} for query {current_query}")
                elif line:
                    self.logger.debug(f"Failed to parse hit line: {line}")

            # Stop hit processing when a new query or blank line is encountered
            if processing_hits and (line.startswith('>>>') or not line):
                processing_hits = False

        # ------------------------------------------------------------------
        # Pass 2: Parse the -m10 alignment section to obtain full descriptions
        # The -m10 query header looks like:
        #   ">>>NM_001001130.3.p1, 518 aa vs protein_db.fasta library"
        # An alignment hit header looks like:
        #   ">>subject_id Full description"
        # ------------------------------------------------------------------
        current_m10_query = None

        for line in lines:
            line = line.strip()

            # -m10 format query header: ">>>query_id, xxx aa vs ..."
            # (no leading digit, unlike the standard "1>>>query_id" header)
            m10_query_match = re.match(r'^>>>([^,<]+),\s+\d+\s+aa\s+vs\s+', line)
            if m10_query_match:
                current_m10_query = m10_query_match.group(1).strip()
                continue

            # Parse the first -m10 hit header while excluding query/control headers (">>>").
            if current_m10_query and line.startswith('>>') and not line.startswith('>>>'):
                hit_header = line[2:].strip()
                header_parts = hit_header.split(maxsplit=1)
                hit_subject_id = header_parts[0] if header_parts else ''
                full_description = header_parts[1].strip() if len(header_parts) == 2 else ''

                # Update the annotation only if it matches the best summary hit for this query.
                if current_m10_query in annotations:
                    best = annotations[current_m10_query]
                    if hit_subject_id == best['subject_id'] and full_description:
                        best['description'] = full_description
                        self.logger.debug(
                            f"Pass 2 - updated description for {current_m10_query}: "
                            f"{full_description[:60]}"
                        )

                # Move on; only the first (best) alignment hit is needed per query.
                current_m10_query = None

        # ------------------------------------------------------------------
        # Write results to output file
        # ------------------------------------------------------------------
        annotation_list = list(annotations.values())

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\tdescription\n")
            for ann in annotation_list:
                f.write(
                    f"{ann['query_id']}\t{ann['subject_id']}\t{ann['identity']}\t"
                    f"{ann['length']}\t{ann['evalue']}\t{ann['bit_score']}\t{ann['description']}\n"
                )

        self.logger.info(f"Parsed {len(annotation_list)} annotations from ggsearch output")
    
    def _parse_ggsearch_line_based(self, content, output_path):
        """Parse ggsearch results line by line (fallback)"""
        annotations = []
        lines = content.split('\n')
        
        current_query = None
        in_hits_section = False
        
        for line in lines:
            line = line.strip()
            
            # Search for the query ID
            if line.startswith('>>>') and ',' in line:
                current_query = line[3:].split(',')[0].strip()
                in_hits_section = False
                continue
            
            # Detect the start of a hit
            if current_query and '>>' in line and 'aa)' in line:
                in_hits_section = True
                
                # Extract hit information
                hit_match = re.search(r'>>(.+?)\s+\((\d+)\s+aa\)', line)
                if hit_match:
                    subject_info = hit_match.group(1).strip()
                    subject_length = hit_match.group(2)
                    
                    # Build a concise annotation summary
                    annotation = {
                        'query_id': current_query,
                        'subject_id': subject_info.split()[0] if subject_info else '',
                        'identity': '',
                        'length': subject_length,
                        'evalue': '',
                        'bit_score': '',
                        'description': subject_info
                    }
                    
                    annotations.append(annotation)
                    # Retrieve only the first hit
                    current_query = None
                    in_hits_section = False
        
        # Save the results to a file
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\tdescription\n")
            for ann in annotations:
                f.write(f"{ann['query_id']}\t{ann['subject_id']}\t{ann['identity']}\t"
                       f"{ann['length']}\t{ann['evalue']}\t{ann['bit_score']}\t{ann['description']}\n")
        
        self.logger.info(f"Line-based parsing found {len(annotations)} annotations")
    
    def _extract_hits(self, query_content):
        """Extract hit information from query content (improved)"""
        hits = []
        
        # Try multiple hit patterns
        hit_patterns = [
            r'>>(.+?)\s+\((\d+)\s+aa\)',  # original
            r'>>(.+?)\s+\((\d+)\s+nt\)',  # nucleotde
            r'>(.+?)\s+Length:\s*(\d+)'    # other formats
        ]
        
        for pattern in hit_patterns:
            matches = list(re.finditer(pattern, query_content))
            if matches:
                self.logger.debug(f"Found {len(matches)} hits with pattern: {pattern}")
                break
        
        if not matches:
            return hits
        
        for match in matches:
            subject_info = match.group(1).strip()
            subject_length = match.group(2)
            
            # Extract the score information
            hit_start = match.start()
            next_match = re.search(r'>>', query_content[hit_start+1:])
            hit_end = next_match.start() + hit_start + 1 if next_match else len(query_content)
            
            hit_section = query_content[hit_start:hit_end]
            
            # Extract the E-value, score, and sequence identity (using multiple patterns)
            evalue = ''
            bit_score = ''
            identity = ''
            
            # Extract the E-value
            evalue_patterns = [
                r'E\(\d+\):\s+([\d.e-]+)',
                r'E-value:\s*([\d.e-]+)',
                r'Expect\s*=\s*([\d.e-]+)'
            ]
            for e_pattern in evalue_patterns:
                e_match = re.search(e_pattern, hit_section)
                if e_match:
                    evalue = e_match.group(1)
                    break
            
            # Extract the score
            score_patterns = [
                r'(\d+\.?\d*)\s+bits',
                r'Score\s*=\s*(\d+\.?\d*)',
                r'bits:\s*(\d+\.?\d*)'
            ]
            for s_pattern in score_patterns:
                s_match = re.search(s_pattern, hit_section)
                if s_match:
                    bit_score = s_match.group(1)
                    break
            
            # Extract the sequence identity
            identity_patterns = [
                r'identity:\s+(\d+\.?\d*)%',
                r'Identity\s*=\s*(\d+\.?\d*)%',
                r'(\d+\.?\d*)%\s+identity'
            ]
            for i_pattern in identity_patterns:
                i_match = re.search(i_pattern, hit_section)
                if i_match:
                    identity = i_match.group(1)
                    break
            
            hit = {
                'subject_id': subject_info.split()[0] if subject_info else '',
                'description': subject_info,
                'length': subject_length,
                'evalue': evalue,
                'bit_score': bit_score,
                'identity': identity
            }
            
            hits.append(hit)
        
        return hits
    
    def classify_genes(self, annotation_path):
        """Classify genes (Multi/Single)"""
        self.logger.info("Start gene classification...")
        
        # Counters for annotation result statistics
        total_annotations = 0
        included_count = 0
        
        # Load the annotation results using subject_id as the key
        subjects = {}
        
        with open(annotation_path, 'r', encoding='utf-8') as f:
            header = f.readline().strip()
            for row_number, line in enumerate(f, 2):
                # Preserve trailing empty fields such as an absent description.
                annotation_line = line.rstrip('\r\n')
                if not annotation_line:
                    continue
                fields = annotation_line.split('\t')
                if len(fields) < 7:
                    if self.protein_transcript_mapping is not None:
                        raise TranscriptMappingError(
                            "Malformed annotation row cannot be reconciled "
                            "with the protein-transcript mapping at "
                            f"{annotation_path}:{row_number}"
                        )
                    continue

                query_id = fields[0]
                subject_id = fields[1]
                description = fields[6]
                if (
                    self.protein_transcript_mapping is not None
                    and (not query_id or not subject_id)
                ):
                    raise TranscriptMappingError(
                        "Annotation row has an empty query_id or subject_id "
                        "and cannot be reconciled with the protein-transcript "
                        f"mapping at {annotation_path}:{row_number}"
                    )
                    
                total_annotations += 1
                    
                # Require every mapped-mode annotation query to remain valid.
                if self.protein_transcript_mapping is not None:
                    if query_id not in self.protein_transcript_mapping:
                        raise TranscriptMappingError(
                            "Annotation query ID is absent from the "
                            "validated protein-transcript mapping: "
                            f"{query_id!r}. Refusing to classify this "
                            "entry as untargetable. Rerun annotation "
                            "with matching query and mapping inputs."
                        )
                    
                # Process for analysis
                included_count += 1
                    
                if subject_id not in subjects:
                    subjects[subject_id] = []
                    
                subjects[subject_id].append({
                    'query_id': query_id,
                    'subject_id': subject_id,
                    'description': description,
                    'line_data': annotation_line,
                })
        
        # Determine Multi/Single status based on subject_id and assign sequential unique_IDs.
        multi_counter = 1
        single_counter = 1
        subject_to_unique_id = {}
        
        for subject_id, entries in subjects.items():
            if len(entries) > 1:
                # Multi: Multiple query_ids match the same subject_id
                unique_id = f"Multi_{multi_counter:010d}"
                subject_to_unique_id[subject_id] = unique_id
                multi_counter += 1
            else:
                # Single: One query_id matches one subject_id
                unique_id = f"Single_{single_counter:010d}"
                subject_to_unique_id[subject_id] = unique_id
                single_counter += 1
        
        # Output the classification results
        classified_path = os.path.join(self.intermediate_dir, 'annotation_classified.txt')
        
        with open(classified_path, 'w', encoding='utf-8') as f:
            f.write("query_id\tsubject_id\tidentity\tlength\tevalue\tbit_score\tdescription\tgene_name\tclassification\tunique_id\n")
            
            # Output the results for each subject_id
            for subject_id, entries in subjects.items():
                unique_id = subject_to_unique_id[subject_id]
                classification = "Multi" if unique_id.startswith("Multi") else "Single"
                
                for entry in entries:
                    # Extract the gene name
                    gene_name = self._extract_gene_name(entry['description'])
                    
                    # Add the classification information and output the results
                    f.write(f"{entry['line_data']}\t{gene_name}\t{classification}\t{unique_id}\n")
        
        # Log the filtering results and summary statistics
        if self.protein_transcript_mapping is not None:
            self.logger.info("=== Mapping-validated annotation results ===")
            self.logger.info(f"Total number of annotations: {total_annotations}")
            self.logger.info(f"Number of entries analyzed: {included_count}")
        else:
            self.logger.info(f"TransDecoder mode: All {total_annotations} annotations will be included in the analysis")
        
        self.logger.info(f"Save the gene classification results: {classified_path}")
        self.logger.info(f"Multi genes: {multi_counter - 1}, Single genes: {single_counter - 1}")
        self.logger.info(f"Total number of genes analyzed: {(multi_counter - 1) + (single_counter - 1)}")
        
        return classified_path
    
    def _extract_gene_name(self, description):
        """Extract gene name from description"""
        # Extract gene names from UniProt descriptions
        # 例: "sp|P12345|GENE_HUMAN Gene name GN=GENE PE=1 SV=1"
        
        # Search for the pattern of "GN="
        gn_match = re.search(r'GN=([^\s]+)', description)
        if gn_match:
            return gn_match.group(1)
        
        # Fallback: use the first word
        parts = description.split()
        if len(parts) > 0:
            return parts[0].split('|')[-1] if '|' in parts[0] else parts[0]
        
        return "unknown"
    
    def extract_gene_sequences(self, classified_path):
        """Extract one locked generation of Multi/Single gene sequences."""
        with annotation_workspace_lock(self.intermediate_dir):
            return self._extract_gene_sequences_locked(classified_path)

    def _extract_gene_sequences_locked(self, classified_path):
        """Extract gene sequences for Multi/Single genes"""
        self.logger.info("Start extracting gene sequences...")
        
        # Load the classified file
        multi_genes = {}
        single_genes = {}
        
        with open(classified_path, 'r', encoding='utf-8') as f:
            header = f.readline().strip()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 10:
                    query_id = fields[0]
                    classification = fields[8]
                    unique_id = fields[9]
                    
                    if classification == "Multi":
                        if unique_id not in multi_genes:
                            multi_genes[unique_id] = []
                        multi_genes[unique_id].append(query_id)
                    else:
                        single_genes[unique_id] = query_id
        
        # Extract sequences from the original transcriptome file
        transcriptome_path = os.path.join(self.intermediate_dir, 'transcriptome.fa')
        transcript_sequences = self._load_transcript_sequences(
            transcriptome_path
        )

        # Validate every Single and Multi member before creating any artifact.
        self._resolve_required_sequences(multi_genes, transcript_sequences)
        self._resolve_required_sequences(single_genes, transcript_sequences)

        # Extract sequences for Multi genes and create separate files
        if multi_genes:
            multi_seq_path = os.path.join(self.intermediate_dir, 'multi_genes.fa')
            self._extract_sequences_by_ids(
                transcriptome_path,
                multi_seq_path,
                multi_genes,
                sequences=transcript_sequences,
            )
            
            # Split multi genes into individual files
            self._create_individual_multi_gene_files(
                transcriptome_path,
                multi_genes,
                sequences=transcript_sequences,
            )
            
            # Execute alignment process
            self._run_alignment_process(multi_genes)
        else:
            clear_multi_gene_artifacts(self.intermediate_dir)
        
        # Extract sequences for Single genes
        if single_genes:
            single_seq_path = os.path.join(self.intermediate_dir, 'single_genes.fa')
            self._extract_sequences_by_ids(
                transcriptome_path,
                single_seq_path,
                single_genes,
                sequences=transcript_sequences,
            )
        else:
            clear_single_gene_artifacts(self.intermediate_dir)
        
        self.logger.info("Gene sequence extraction completed successfully")
        
        return {
            'multi_genes': multi_genes,
            'single_genes': single_genes
        }
    
    @staticmethod
    def _load_transcript_sequences(fasta_path):
        """Load transcript sequences keyed by the first FASTA header token."""
        return read_fasta_sequences(fasta_path, 'transcript')

    def _extract_sequences_by_ids(
        self,
        fasta_path,
        output_path,
        gene_dict,
        sequences=None,
    ):
        """Extract sequences for specified IDs"""
        if sequences is None:
            sequences = self._load_transcript_sequences(fasta_path)

        resolved_sequences = self._resolve_required_sequences(
            gene_dict, sequences
        )

        # Output the sequence for the specified ID only after every ID resolves.
        with open(output_path, 'w') as f:
            for unique_id, query_ids in gene_dict.items():
                if isinstance(query_ids, list):
                    # For Multi gene
                    for query_id, sequence_id in resolved_sequences[unique_id]:
                        f.write(
                            f">{unique_id}_{query_id}\n"
                            f"{sequences[sequence_id]}\n"
                        )
                else:
                    # For Single gene
                    _, sequence_id = resolved_sequences[unique_id][0]
                    f.write(f">{unique_id}\n{sequences[sequence_id]}\n")

    def _resolve_required_sequences(self, gene_dict, sequences):
        """Resolve every classified query before publishing sequence FASTAs."""
        resolved = {}
        missing = []

        for unique_id, query_ids in gene_dict.items():
            query_id_list = (
                query_ids if isinstance(query_ids, list) else [query_ids]
            )
            resolved[unique_id] = []
            for query_id in query_id_list:
                sequence_id = self._resolve_sequence_id(query_id)
                if sequence_id not in sequences:
                    missing.append(f"{query_id}->{sequence_id}")
                    continue
                resolved[unique_id].append((query_id, sequence_id))

        if missing:
            examples = ', '.join(missing[:5])
            if len(missing) > 5:
                examples += f", ... (+{len(missing) - 5})"
            if self.protein_transcript_mapping is not None:
                corrective_action = (
                    "correct the protein-transcript mapping and rerun "
                    "annotation"
                )
            else:
                corrective_action = (
                    "verify that query.pep uses TransDecoder-style "
                    "<transcript_id>.p<number> identifiers matching "
                    "transcriptome.fa, then rerun gene prediction and annotation"
                )
            raise TranscriptMappingError(
                "Classified protein query ID(s) could not be resolved to "
                "transcriptome FASTA records: "
                f"{examples}. Refusing to create incomplete gene FASTAs; "
                f"{corrective_action}."
            )

        return resolved
    
    def _resolve_sequence_id(self, query_id):
        """
        Resolve actual sequence ID from query ID
        
        1. Use the protein-to-transcript mapping if available
        2. When no mapping is available, use the conventional approach (e.g., remove suffixes such as .p1)
        
        Args:
            query_id (str): Query ID (typically a protein ID)
            
        Returns:
            str: The ID actually used in the FASTA file, or None
        """
        # If a mapping file has been loaded
        if self.protein_transcript_mapping is not None:
            # protein_id -> transcript_id Attempt the mapping
            if query_id in self.protein_transcript_mapping:
                transcript_id = self.protein_transcript_mapping[query_id]
                self.logger.debug(f"Apply the mapping: {query_id} -> {transcript_id}")
                return transcript_id
            raise TranscriptMappingError(
                "Protein query ID is absent from the loaded "
                f"protein-transcript mapping: {query_id!r}"
            )
        
        # Conventional processing when no mapping is available
        # Remove suffixes such as .p1 from query_id and use the result as the transcript ID
        base_id = re.sub(r'\.p\d+$', '', query_id)
        self.logger.debug(f"従来の処理適用: {query_id} -> {base_id}")
        return base_id
    
    def _validate_existing_protein_mode_requirements(self):
        """
        Validate input requirements for existing protein mode
        
        When using existing protein sequences, ensure that a mapping file is provided
        """
        input_mode = determine_protein_input_mode(self.config)

        if input_mode == EXISTING_PROTEIN_MODE:
            self.logger.info("Existing protein mode detected")
            
            # Check whether the mapping file exists
            mapping_path = self.config['input'].get('protein_transcript_mapping', '')
            
            if not mapping_path or not mapping_path.strip():
                raise ValueError(
                    "In existing protein mode, input.protein_transcript_mapping "
                    "is required. Specify a headerless "
                    "transcript_id<TAB>protein_id mapping file."
                )
            
            if not os.path.exists(mapping_path):
                raise FileNotFoundError(
                    f"The specified mapping file was not found: {mapping_path}\n"
                    "Existing protein mode requires a readable "
                    "input.protein_transcript_mapping file."
                )
            
        else:
            self.logger.info("Use TransDecoder mode (no mapping file required)")

        return input_mode
    
    def _create_individual_multi_gene_files(
        self,
        transcriptome_path,
        multi_genes,
        sequences=None,
    ):
        """Split multi genes into individual files"""
        self.logger.info("Split multi genes into individual files中...")
        
        # Create the multi_genes directory
        multi_genes_dir = os.path.join(self.intermediate_dir, 'multi_genes')
        multi_genes_dir_path = Path(multi_genes_dir)
        multi_genes_dir_path.mkdir(parents=True, exist_ok=True)

        # Do not let groups from an older classification enter this run.
        expected_group_ids = set(multi_genes)
        for previous_group in multi_genes_dir_path.glob('*.fa'):
            if previous_group.stem not in expected_group_ids:
                previous_group.unlink()
        
        if sequences is None:
            sequences = self._load_transcript_sequences(transcriptome_path)

        resolved_sequences = self._resolve_required_sequences(
            multi_genes, sequences
        )

        # Create a separate file for each Multi gene
        for unique_id, query_ids in multi_genes.items():
            output_file = os.path.join(multi_genes_dir, f"{unique_id}.fa")
            
            with open(output_file, 'w') as f:
                for query_id, sequence_id in resolved_sequences[unique_id]:
                    f.write(
                        f">{unique_id}_{query_id}\n"
                        f"{sequences[sequence_id]}\n"
                    )
            
            self.logger.debug(f"Create individual files for Multi genes: {output_file}")
        
        self.logger.info(f"Individual files for Multi genes created successfully: {len(multi_genes)} file")
    
    def _run_alignment_process(self, expected_multi_ids=None):
        """Run cached, bounded-parallel alignment for all Multi genes."""
        runner = MultiGeneAlignmentRunner(
            self.intermediate_dir,
            self.config,
            self.logger,
            run_command=subprocess.run,
        )
        return runner.run(expected_multi_ids)
    
    def run(self):
        """Execute one serialized annotation generation."""
        with annotation_workspace_lock(self.intermediate_dir):
            return self._run_locked()

    def _run_locked(self):
        """Execute entire annotation processing"""
        try:
            # Check the results from the previous step
            query_path = os.path.join(self.intermediate_dir, 'query.pep')
            if not os.path.exists(query_path):
                raise FileNotFoundError(f"Predicted protein sequence file not found: {query_path}")
            
            # Input validation: Check the mapping file requirement for existing protein mode
            input_mode = self._validate_existing_protein_mode_requirements()

            transcriptome_path = os.path.join(
                self.intermediate_dir, 'transcriptome.fa'
            )
            if input_mode == EXISTING_PROTEIN_MODE:
                # Validate the exact staged inputs before the expensive search.
                self.load_protein_transcript_mapping()
                self.validate_protein_transcript_mapping(
                    query_path,
                    transcriptome_path,
                )
            else:
                # A configured mapping belongs only to existing protein mode.
                self.protein_transcript_mapping = None
            
            # Prepare protein database for GO analysis
            db_path = self.prepare_protein_db_for_goanalysis()
            
            # FASTA search
            search_results = self.run_fasta_search(query_path, db_path)
            
            # Parsing the results
            annotation_table = self.parse_fasta_results(search_results)
            
            # Gene classification
            classified_table = self.classify_genes(annotation_table)
            
            # Gene sequence extraction
            gene_sequences = self.extract_gene_sequences(classified_table)
            
            # Save the results
            results = {
                'db_path': db_path,
                'search_results': search_results,
                'annotation_table': annotation_table,
                'classified_table': classified_table,
                'gene_sequences': gene_sequences
            }
            
            # Save the results to a file
            results_file = os.path.join(self.intermediate_dir, 'annotation_results.txt')
            with open(results_file, 'w', encoding='utf-8') as f:
                f.write("=== Annotation results ===\n\n")
                f.write(f"UniProt database: {db_path}\n")
                f.write(f"FASTA search result: {search_results}\n")
                f.write(f"Annotation table: {annotation_table}\n")
                f.write(f"Classification table: {classified_table}\n")
                f.write(f"Multi genes: {len(gene_sequences['multi_genes'])}\n")
                f.write(f"Single genes: {len(gene_sequences['single_genes'])}\n")
            
            self.logger.info(f"Save the annotation results: {results_file}")
            self.logger.info("Annotation completed successfully")
            
            return results
            
        except Exception as e:
            self.logger.error(f"An error occurred during annotation: {e}")
            raise
