#!/usr/bin/env python3
"""
Gene Prediction Module
ORF prediction using TransDecoder and quality assessment using BUSCO
"""

import os
import subprocess
import logging
from pathlib import Path
import gzip
import shutil

from .input_mode import EXISTING_PROTEIN_MODE, determine_protein_input_mode


class GenePrediction:
    """Gene prediction processing class"""

    TRANSDECODER_DEFAULT_MIN_PROTEIN_LEN = 100
    
    def __init__(self, config):
        """
        Initialize the gene prediction module
        
        Args:
            config (dict): Configuration dictionary
        """
        self.config = config
        self.logger = logging.getLogger(__name__)
        self.intermediate_dir = self.config['output']['intermediate']
        
    def determine_input_mode(self):
        """Determine input mode (TransDecoder or existing protein sequences)"""
        transcriptome_path = self.config['input'].get('transcriptome', '')
        
        # Transcriptome is always required (tWhEEL is a transcriptome-based pipeline)
        if not transcriptome_path or not os.path.exists(transcriptome_path):
            raise FileNotFoundError("Transcriptome file is required for tWhEEL pipeline")

        input_mode = determine_protein_input_mode(self.config)
        if input_mode == EXISTING_PROTEIN_MODE:
            self.logger.info("Existing protein sequences mode - skipping TransDecoder")
            self.logger.info("Transcriptome data will be used in downstream analyses")
        else:
            self.logger.info("Standard TransDecoder mode - ORF prediction from transcriptome")
        return input_mode
    
    def prepare_input(self):
        """Prepare input transcriptome file"""
        transcriptome_path = self.config['input']['transcriptome']
        
        # Output file path
        output_path = os.path.join(self.intermediate_dir, 'transcriptome.fa')
        
        # Decompress and copy file
        if transcriptome_path.endswith('.gz'):
            self.logger.info("Decompressing transcriptome file...")
            with gzip.open(transcriptome_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            self.logger.info("Copying transcriptome file...")
            shutil.copy2(transcriptome_path, output_path)
        
        # Format sequences to one-line format
        self._format_fasta_oneline(output_path)
        
        self.logger.info(f"Input file preparation complete: {output_path}")
        return output_path
    
    def prepare_existing_protein_sequences(self):
        """Prepare existing protein sequence file"""
        protein_path = self.config['input']['protein_sequences']
        
        # Output file path
        output_path = os.path.join(self.intermediate_dir, 'query.pep')
        
        self.logger.info("Preparing existing protein sequence file...")
        
        # Decompress and copy file
        if protein_path.endswith('.gz'):
            self.logger.info("Decompressing protein sequence file...")
            with gzip.open(protein_path, 'rb') as f_in:
                with open(output_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            self.logger.info("Copying protein sequence file...")
            shutil.copy2(protein_path, output_path)
        
        # Format sequences to one-line format
        self._format_fasta_oneline(output_path)
        
        self.logger.info(f"Protein sequence file preparation complete: {output_path}")
        return output_path
    
    def _format_fasta_oneline(self, fasta_path):
        """Convert FASTA file to one-line format"""
        temp_path = fasta_path + '.temp'
        
        with open(fasta_path, 'r') as f_in:
            with open(temp_path, 'w') as f_out:
                sequence = ""
                header = ""
                
                for line in f_in:
                    line = line.strip()
                    if line.startswith('>'):
                        if header and sequence:
                            f_out.write(f"{header}\n{sequence}\n")
                        header = line
                        sequence = ""
                    else:
                        sequence += line
                
                # Output last sequence
                if header and sequence:
                    f_out.write(f"{header}\n{sequence}\n")
        
        # Replace original file
        shutil.move(temp_path, fasta_path)
        self.logger.debug("FASTA file converted to one-line format")

    def _get_transdecoder_min_protein_len(self):
        """Return a validated TransDecoder minimum protein length."""
        value = self.config.get('transdecoder', {}).get(
            'min_protein_len', self.TRANSDECODER_DEFAULT_MIN_PROTEIN_LEN
        )
        if isinstance(value, bool):
            raise ValueError(
                'transdecoder.min_protein_len must be a positive integer'
            )
        if isinstance(value, int):
            parsed_value = value
        elif isinstance(value, str):
            try:
                parsed_value = int(value.strip())
            except ValueError as exc:
                raise ValueError(
                    'transdecoder.min_protein_len must be a positive integer'
                ) from exc
        else:
            raise ValueError(
                'transdecoder.min_protein_len must be a positive integer'
            )
        if parsed_value <= 0:
            raise ValueError(
                'transdecoder.min_protein_len must be a positive integer'
            )
        return parsed_value
    
    def run_transdecoder(self, input_path):
        """ORF prediction using TransDecoder"""
        self.logger.info("Begin ORF prediction using TransDecoder...")
        min_protein_len = self._get_transdecoder_min_protein_len()
        self.logger.info(
            f"TransDecoder minimum protein length: {min_protein_len} aa"
        )
        
        # Create a dedicated working directory for TransDecoder
        transdecoder_workdir = os.path.join(self.intermediate_dir, 'transdecoder_output')
        Path(transdecoder_workdir).mkdir(parents=True, exist_ok=True)
        
        # Copy the input files to the working directory
        input_filename = os.path.basename(input_path)
        work_input_path = os.path.join(transdecoder_workdir, input_filename)
        shutil.copy2(input_path, work_input_path)
        
        # Docker settings
        docker_image = self.config['transdecoder']['docker_image']
        
        # Run TransDecoder.LongOrfs from within the working directory
        cmd_longorfs = [
            "docker", "run", 
            "-v", f"{os.path.abspath(transdecoder_workdir)}:/work",
            "-w", "/work",
            "--rm",
            docker_image,
            "TransDecoder.LongOrfs",
            "-t", input_filename,
            "-m", str(min_protein_len),
        ]
        
        self.logger.debug(f"TransDecoder.LongOrfs実行コマンド: {' '.join(cmd_longorfs)}")
        result = subprocess.run(cmd_longorfs, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"TransDecoder.LongOrfs実行エラー: {result.stderr}")
            raise RuntimeError("TransDecoder.LongOrfs実行に失敗しました")
        
        # Run TransDecoder.Predict from within the working directory
        cmd_predict = [
            "docker", "run",
            "-v", f"{os.path.abspath(transdecoder_workdir)}:/work",
            "-w", "/work", 
            "--rm",
            docker_image,
            "TransDecoder.Predict",
            "-t", input_filename
        ]
        
        self.logger.debug(f"TransDecoder.LongOrfs command: {' '.join(cmd_predict)}")
        result = subprocess.run(cmd_predict, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"TransDecoder.Predict execution error: {result.stderr}")
            raise RuntimeError("TransDecoder.Predict failed to run")
        
        # Copy the predicted protein sequence file to the final output directory
        protein_filename = f"{input_filename}.transdecoder.pep"
        source_protein_path = os.path.join(transdecoder_workdir, protein_filename)
        query_path = os.path.join(self.intermediate_dir, 'query.pep')
        
        if os.path.exists(source_protein_path):
            shutil.copy2(source_protein_path, query_path)
            self.logger.info(f"Save the predicted protein sequences: {query_path}")
            self.logger.info(f"TransDecoder working directory: {transdecoder_workdir}")
        else:
            raise FileNotFoundError(f"TransDecoder output file not found: {source_protein_path}")
        
        return query_path
    
    def get_sequence_stats(self, protein_path):
        """Retrieve statistics for the protein sequences"""
        self.logger.info("Retrieving protein sequence statistics...")
        
        cmd = [
            "docker", "run",
            "-v", f"{os.getcwd()}:{os.getcwd()}",
            "-w", os.getcwd(),
            "--rm",
            "quay.io/biocontainers/seqkit:2.8.0--h9ee0642_0",
            "seqkit", "stats", protein_path
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.warning(f"Error retrieving statistics with seqkit: {result.stderr}")
            return None
        
        self.logger.info("Protein sequence statistics:")
        self.logger.info(result.stdout)
        
        return result.stdout
    
    def run_busco(self, protein_path):
        """Run the BUSCO assessment"""
        self.logger.info("Start the BUSCO assessment...")
        
        # BUSCO output directory name (name only, not the full path)
        busco_output_name = 'busco_output'
        busco_output_path = os.path.join(self.intermediate_dir, busco_output_name)
        
        # Get the relative path to protein_path.
        protein_filename = os.path.basename(protein_path)
        
        cmd = [
            "docker", "run",
            "--platform", "linux/amd64",  # for Apple Silicon
            "-v", f"{os.path.abspath(self.intermediate_dir)}:/busco_work",
            "-w", "/busco_work",
            "--rm",
            self.config['busco']['docker_image'],
            "busco",
            "-m", "protein",
            "-i", protein_filename,  # Relative path
            "-o", busco_output_name,  # only directory name
            "-l", self.config['busco']['lineage'],
            "-c", str(self.config['busco']['cpu'])
        ]
        
        self.logger.debug(f"BUSCO command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            self.logger.error(f"BUSCO execution error: {result.stderr}")
            raise RuntimeError("BUSCO assessment failed")
        
        # Reading BUSCO results
        busco_summary = self._parse_busco_results(busco_output_path)
        
        self.logger.info("BUSCO assessment results:")
        self.logger.info(busco_summary)
        
        return busco_summary
    
    def _parse_busco_results(self, busco_output):
        """Analysis of BUSCO results"""
        summary_files = list(Path(busco_output).glob("**/short_summary*.txt"))
        
        if not summary_files:
            self.logger.warning("BUSCO short summary file not found")
            return "Failed to analyze the BUSCO results"
        
        summary_file = summary_files[0]
        
        try:
            with open(summary_file, 'r') as f:
                content = f.read()
            return content
        except Exception as e:
            self.logger.error(f"Error reading BUSCO results: {e}")
            return "Unable to read the BUSCO results"
    
    def run(self):
        """Run the complete gene prediction workflow"""
        try:
            # Determine the input mode
            input_mode = self.determine_input_mode()
            
            # Prepare the transcriptome file in all modes (a core requirement of tWhEEL)
            transcriptome_path = self.prepare_input()
            
            if input_mode == "existing_protein":
                # Use the existing protein sequences
                protein_path = self.prepare_existing_protein_sequences()
                input_path = transcriptome_path  # transcriptomeパスを記録
                transdecoder_skipped = True
                self.logger.info("Use the existing protein sequences and skip TransDecoder")
                self.logger.info("The transcript data will be used in downstream analyses, such as guide RNA design")
                
            elif input_mode == "transdecoder":
                # Perform ORF prediction using TransDecoder
                input_path = transcriptome_path
                protein_path = self.run_transdecoder(input_path)
                transdecoder_skipped = False
                
            else:
                raise ValueError(f"Unknown input mode: {input_mode}")
            
            # Retrieve sequence statistics
            stats = self.get_sequence_stats(protein_path)
            
            # Perform BUSCO assessment
            busco_result = self.run_busco(protein_path)
            
            # Save the results
            results = {
                'input_mode': input_mode,
                'input_path': input_path,
                'protein_path': protein_path,
                'transdecoder_skipped': transdecoder_skipped,
                'sequence_stats': stats,
                'busco_result': busco_result
            }
            
            # Save the results to a file
            results_file = os.path.join(self.intermediate_dir, 'gene_prediction_results.txt')
            with open(results_file, 'w', encoding='utf-8') as f:
                f.write("=== Gene prediction results ===\n\n")
                f.write(f"Input mode: {input_mode}\n")
                f.write(f"Input file: {input_path}\n")
                f.write(f"Protein sequences: {protein_path}\n")
                if transdecoder_skipped:
                    f.write("TransDecoder: Skipped (using existing protein sequences)\n")
                else:
                    f.write("TransDecoder: Completed\n")
                f.write("\n=== Sequence statistics ===\n")
                f.write(stats if stats else "Failed to retrieve sequence statistics\n")
                f.write("\n=== BUSCO assessment results ===\n")
                f.write(busco_result)
            
            self.logger.info(f"Save the gene prediction results: {results_file}")
            self.logger.info("Gene prediction completed successfully")
            
            return results
            
        except Exception as e:
            self.logger.error(f"An error occurred during gene prediction: {e}")
            raise
