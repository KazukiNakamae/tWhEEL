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
        
        mapping_dict = {}
        try:
            with open(mapping_path, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):  # 空行やコメント行をスキップ
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        transcript_id = parts[0].strip()
                        protein_id = parts[1].strip()
                        # Map protein_id to transcript_id
                        mapping_dict[protein_id] = transcript_id
                    else:
                        self.logger.warning(f"Mapping file entries {line_num}: Invalid format - {line}")
            
            self.logger.info(f"Build the mapping dictionary: {len(mapping_dict)} entry")
            self.protein_transcript_mapping = mapping_dict
            return mapping_dict
            
        except Exception as e:
            self.logger.error(f"Error reading the mapping file: {e}")
            raise RuntimeError(f"Failed to read the mapping file: {mapping_path}")
    
    def run_fasta_search(self, query_path, db_path):
        """Execute FASTA search"""
        # Output file path
        output_path = os.path.join(self.intermediate_dir, 'ggsearch_results.txt')
        
        # Check for existing files
        if os.path.exists(output_path):
            # Check the file size to ensure the file is not empty
            if os.path.getsize(output_path) > 0:
                self.logger.info(f"Existing FASTA search results found: {output_path}")
                self.logger.info("Skipping the FASTA search")
                return output_path
            else:
                self.logger.warning(f"The existing FASTA search results file is empty. Rerunning the search: {output_path}")
        
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
        
        # 標準出力をファイルに保存
        with open(output_path, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"FASTA search execution error: {result.stderr}")
            raise RuntimeError("FASTA search failed")
        
        self.logger.info(f"Save the FASTA search results: {output_path}")
        return output_path
    
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
        2. Machine-readable (-m10) alignment section with ">>sp|..." headers - descriptions are FULL
           (e.g., "Zinc finger protein 708 OS=Homo sapiens OX=9606 GN=ZNF708 PE=1 SV=6")

        This method performs a two-pass parse:
        - Pass 1: Reads the standard "best scores" section to obtain subject_id, length, bits, evalue
          for the best (first) hit per query.
        - Pass 2: Reads the -m10 alignment header lines (">>sp|...") to obtain the FULL description
          including the GN= field, then updates the annotations collected in Pass 1.
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
        # The best-score hit line looks like:
        #   "sp|P17019|ZN708_HUMAN Zinc finger protein 708 OS=H ( 563) 1854 97.9 0"
        # ------------------------------------------------------------------
        current_query = None
        processing_hits = False

        for line in lines:
            line = line.strip()

            # Standard-format query header: "\d+>>>query_id ... - xxx aa"
            query_match = re.match(r'^\d+>>>(.+?)\s+.+?\s+-\s+\d+\s+aa', line)
            if query_match:
                current_query = query_match.group(1).strip()
                processing_hits = False
                self.logger.debug(f"Found query: {current_query}")
                continue

            # "The best scores are:" marks the start of the hit summary for this query
            if "The best scores are:" in line:
                processing_hits = True
                continue

            # First hit line in the "best scores" section (truncated description)
            if processing_hits and current_query and current_query not in annotations and line.startswith('sp|'):
                hit_match = re.match(
                    r'sp\|([^|]+)\|([^\s]+)\s+(.+?)\s+\(\s*(\d+)\)\s+(\d+)\s+([\d.]+)\s+([\d.e-]+)',
                    line
                )
                if hit_match:
                    subject_acc  = hit_match.group(1)
                    subject_name = hit_match.group(2)
                    description  = hit_match.group(3).strip()   # Truncated; will be replaced in Pass 2
                    length       = hit_match.group(4)
                    bits         = hit_match.group(6)
                    evalue       = hit_match.group(7)

                    annotations[current_query] = {
                        'query_id':   current_query,
                        'subject_id': f"sp|{subject_acc}|{subject_name}",
                        'subject_acc': subject_acc,
                        'identity':   '',
                        'length':     length,
                        'evalue':     evalue,
                        'bit_score':  bits,
                        'description': description,  # Truncated (fallback if Pass 2 fails)
                    }
                    processing_hits = False
                    self.logger.debug(f"Pass 1 - best hit: {subject_acc} for query {current_query}")
                else:
                    self.logger.debug(f"Failed to parse hit line: {line}")

            # Stop hit processing when a new query or blank line is encountered
            if processing_hits and (line.startswith('>>>') or not line):
                processing_hits = False

        # ------------------------------------------------------------------
        # Pass 2: Parse the -m10 alignment section to obtain full descriptions
        # The -m10 query header looks like:
        #   ">>>NM_001001130.3.p1, 518 aa vs protein_db.fasta library"
        # The alignment hit header looks like:
        #   ">>sp|P17019|ZN708_HUMAN Zinc finger protein 708 OS=Homo sapiens OX=9606 GN=ZNF708 PE=1 SV=6"
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

            # -m10 format alignment hit header: ">>sp|ACC|NAME Full description"
            if current_m10_query and line.startswith('>>sp|'):
                align_match = re.match(r'^>>sp\|([^|]+)\|([^\s]+)\s+(.+)$', line)
                if align_match:
                    hit_acc          = align_match.group(1)
                    full_description = align_match.group(3).strip()

                    # Update the annotation only if it matches our best hit for this query
                    if current_m10_query in annotations:
                        best = annotations[current_m10_query]
                        if hit_acc == best['subject_acc']:
                            best['description'] = full_description
                            self.logger.debug(
                                f"Pass 2 - updated description for {current_m10_query}: "
                                f"{full_description[:60]}"
                            )

                    # Move on; we only need the first (best) hit per query
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
        excluded_count = 0
        included_count = 0
        
        # Load the annotation results using subject_id as the key
        subjects = {}
        
        with open(annotation_path, 'r', encoding='utf-8') as f:
            header = f.readline().strip()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 7:
                    query_id = fields[0]
                    subject_id = fields[1]
                    description = fields[6]
                    
                    total_annotations += 1
                    
                    # Filtering when the mapping file has been loaded in existing protein mode
                    if self.protein_transcript_mapping is not None:
                        if query_id not in self.protein_transcript_mapping:
                            # Exclude IDs that are not present in the mapping file
                            excluded_count += 1
                            self.logger.debug(f"Excluded ID (not found in the mapping file): {query_id}")
                            continue
                    
                    # Process for analysis
                    included_count += 1
                    
                    if subject_id not in subjects:
                        subjects[subject_id] = []
                    
                    subjects[subject_id].append({
                        'query_id': query_id,
                        'subject_id': subject_id,
                        'description': description,
                        'line_data': line.strip()
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
            self.logger.info(f"=== Mapping-based filtering results ===")
            self.logger.info(f"Total number of annotations: {total_annotations}")
            self.logger.info(f"Number of entries analyzed: {included_count}")
            self.logger.info(f"Number of excluded entries: {excluded_count}")
            if total_annotations > 0:
                inclusion_rate = (included_count / total_annotations) * 100
                self.logger.info(f"Analysis inclusion rate: {inclusion_rate:.2f}%")
            
            if excluded_count > 0:
                self.logger.warning(f"{excluded_count} ID was excluded because it was not found in the mapping file")
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
        
        # Extract sequences for Multi genes and create separate files
        if multi_genes:
            multi_seq_path = os.path.join(self.intermediate_dir, 'multi_genes.fa')
            self._extract_sequences_by_ids(transcriptome_path, multi_seq_path, multi_genes)
            
            # Split multi genes into individual files
            self._create_individual_multi_gene_files(transcriptome_path, multi_genes)
            
            # Execute alignment process
            self._run_alignment_process()
        
        # Extract sequences for Single genes
        if single_genes:
            single_seq_path = os.path.join(self.intermediate_dir, 'single_genes.fa')
            self._extract_sequences_by_ids(transcriptome_path, single_seq_path, single_genes)
        
        self.logger.info("Gene sequence extraction completed successfully")
        
        return {
            'multi_genes': multi_genes,
            'single_genes': single_genes
        }
    
    def _extract_sequences_by_ids(self, fasta_path, output_path, gene_dict):
        """Extract sequences for specified IDs"""
        # Load all sequences
        sequences = {}
        
        with open(fasta_path, 'r') as f:
            current_id = None
            current_seq = ""
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        sequences[current_id] = current_seq
                    current_id = line[1:].split()[0]  # Retrieve only the first part
                    current_seq = ""
                else:
                    current_seq += line
            
            # Process the last sequence
            if current_id and current_seq:
                sequences[current_id] = current_seq
        
        # Output the sequence for the specified ID
        with open(output_path, 'w') as f:
            for unique_id, query_ids in gene_dict.items():
                if isinstance(query_ids, list):
                    # For Multi gene
                    for query_id in query_ids:
                        sequence_id = self._resolve_sequence_id(query_id)
                        if sequence_id and sequence_id in sequences:
                            f.write(f">{unique_id}_{query_id}\n{sequences[sequence_id]}\n")
                        else:
                            self.logger.warning(f"Sequence not found (Multi): {query_id} -> {sequence_id}")
                else:
                    # For Single gene
                    query_id = query_ids
                    sequence_id = self._resolve_sequence_id(query_id)
                    if sequence_id and sequence_id in sequences:
                        f.write(f">{unique_id}\n{sequences[sequence_id]}\n")
                    else:
                        self.logger.warning(f"Sequence not found (Single): {query_id} -> {sequence_id}")
    
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
        if self.protein_transcript_mapping:
            # protein_id -> transcript_id Attempt the mapping
            if query_id in self.protein_transcript_mapping:
                transcript_id = self.protein_transcript_mapping[query_id]
                self.logger.debug(f"Apply the mapping: {query_id} -> {transcript_id}")
                return transcript_id
            else:
                self.logger.debug(f"Mapping not found: {query_id}")
        
        # Conventional processing when no mapping is available
        # Remove suffixes such as .p1 from query_id and use the result as the transcript ID
        base_id = query_id.split('.p')[0] if '.p' in query_id else query_id
        self.logger.debug(f"従来の処理適用: {query_id} -> {base_id}")
        return base_id
    
    def _validate_existing_protein_mode_requirements(self):
        """
        Validate input requirements for existing protein mode
        
        When using existing protein sequences, ensure that a mapping file is provided
        """
        # Determine whether the current mode is existing protein mode
        protein_sequences = self.config['input'].get('protein_sequences', '')
        is_existing_protein_mode = bool(protein_sequences and protein_sequences.strip())
        
        if is_existing_protein_mode:
            self.logger.info("Existing protein mode detected")
            
            # Check whether the mapping file exists
            mapping_path = self.config['input'].get('protein_transcript_mapping', '')
            
            if not mapping_path or not mapping_path.strip():
                raise ValueError(
                    "In existing protein mode, the protein_transcript_mapping parameter is required"
                    "Please specify the protein ID–to–transcript ID mapping file"
                )
            
            if not os.path.exists(mapping_path):
                raise FileNotFoundError(
                    f"The specified mapping file was not found: {mapping_path}\n"
                    f"In existing protein mode, a protein_transcript_mapping file is required"
                )
            
            # Basic validation of the file format
            try:
                with open(mapping_path, 'r', encoding='utf-8') as f:
                    lines = f.readlines()[:10]  # Check the first 10 lines
                    
                valid_lines = 0
                for line_num, line in enumerate(lines, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        valid_lines += 1
                    else:
                        self.logger.warning(f"Mapping file line {line_num}: The file may not be in TSV format - {line}")
                
                if valid_lines == 0:
                    raise ValueError(
                        f"No valid data found in the mapping file: {mapping_path}\n"
                        f"Please ensure that the file is tab-delimited (TSV format)。"
                    )
                
                self.logger.info(f"Mapping file format validation completed: {valid_lines} Valid entries detected")
                    
            except Exception as e:
                raise RuntimeError(f"Failed to load or validate the mapping file: {mapping_path}\nError: {e}")
            
        else:
            self.logger.info("Use TransDecoder mode (no mapping file required)")
    
    def _create_individual_multi_gene_files(self, transcriptome_path, multi_genes):
        """Split multi genes into individual files"""
        self.logger.info("Split multi genes into individual files中...")
        
        # Create the multi_genes directory
        multi_genes_dir = os.path.join(self.intermediate_dir, 'multi_genes')
        Path(multi_genes_dir).mkdir(parents=True, exist_ok=True)
        
        # Load all sequences
        sequences = {}
        with open(transcriptome_path, 'r') as f:
            current_id = None
            current_seq = ""
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id and current_seq:
                        sequences[current_id] = current_seq
                    current_id = line[1:].split()[0]  # Retrieve only the first part
                    current_seq = ""
                else:
                    current_seq += line
            
            # Process the last sequence
            if current_id and current_seq:
                sequences[current_id] = current_seq
        
        # Create a separate file for each Multi gene
        for unique_id, query_ids in multi_genes.items():
            output_file = os.path.join(multi_genes_dir, f"{unique_id}.fa")
            
            with open(output_file, 'w') as f:
                for query_id in query_ids:
                    sequence_id = self._resolve_sequence_id(query_id)
                    if sequence_id and sequence_id in sequences:
                        f.write(f">{unique_id}_{query_id}\n{sequences[sequence_id]}\n")
                    else:
                        self.logger.warning(f"Sequence not found (individual Multi gene): {query_id} -> {sequence_id}")
            
            self.logger.debug(f"Create individual files for Multi genes: {output_file}")
        
        self.logger.info(f"Individual files for Multi genes created successfully: {len(multi_genes)} file")
    
    def _run_alignment_process(self):
        """Execute alignment process"""
        self.logger.info("Start Multi gene alignment...")
        
        multi_genes_dir = os.path.join(self.intermediate_dir, 'multi_genes')
        multi_genes_aligned_dir = os.path.join(self.intermediate_dir, 'multi_genes_aligned')
        
        # Create the output directory
        Path(multi_genes_aligned_dir).mkdir(parents=True, exist_ok=True)
        
        # Building the Docker image
        docker_build_dir = "src/docker/run_alignment_and_mask"
        self.logger.info("Building run_alignment_and_maskのDocker image...")
        
        build_cmd = [
            "sudo", "docker", "build",
            "-t", "run_alignment_and_mask/v1.0",
            docker_build_dir
        ]
        
        build_result = subprocess.run(build_cmd, capture_output=True, text=True)
        if build_result.returncode != 0:
            self.logger.error(f"Docker build error: {build_result.stderr}")
            raise RuntimeError("Failed to build the Docker image")
        
        self.logger.info("Docker image built successfully")
        
        # Execute alignment process（Specify the path relative to intermediate_dir）
        multi_genes_rel_dir = os.path.relpath(multi_genes_dir, self.intermediate_dir)
        multi_genes_aligned_rel_dir = os.path.relpath(multi_genes_aligned_dir, self.intermediate_dir)
        
        alignment_cmd = [
            "docker", "run",
            "-m", "4G",
            "-v", f"{os.path.abspath(self.intermediate_dir)}:/data",
            "-w", "/data",
            "--rm",
            "run_alignment_and_mask/v1.0",
            "python", "/app/run_alignment_and_mask.py",
            "-i", multi_genes_rel_dir,
            "-o", multi_genes_aligned_rel_dir
        ]
        
        self.logger.debug(f"Alignment command: {' '.join(alignment_cmd)}")
        
        alignment_result = subprocess.run(alignment_cmd, capture_output=True, text=True)
        if alignment_result.returncode != 0:
            self.logger.error(f"Alignment error: {alignment_result.stderr}")
            raise RuntimeError("Alignment failed")
        
        # Log the execution results
        if alignment_result.stdout:
            self.logger.info(f"Alignment output:\n{alignment_result.stdout}")
        
        self.logger.info("Multi gene alignment completed successfully")
        
        # Check the number of output files
        aligned_files = list(Path(multi_genes_aligned_dir).glob("*_aligned_consensus.fa"))
        self.logger.info(f"Number of aligned files generated: {len(aligned_files)}")
    
    def run(self):
        """Execute entire annotation processing"""
        try:
            # Check the results from the previous step
            query_path = os.path.join(self.intermediate_dir, 'query.pep')
            if not os.path.exists(query_path):
                raise FileNotFoundError(f"Predicted protein sequence file not found: {query_path}")
            
            # Input validation: Check the mapping file requirement for existing protein mode
            self._validate_existing_protein_mode_requirements()
            
            # Load the protein-to-transcript mapping file
            self.load_protein_transcript_mapping()
            
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
