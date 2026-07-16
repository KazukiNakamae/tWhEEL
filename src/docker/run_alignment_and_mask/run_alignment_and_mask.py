#!/usr/bin/env python3
"""
Align Multi-gene transcript sequences and generate consensus sequences.
"""

import os
import sys
import argparse
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentProcessor:
    """Align transcript groups and generate one consensus per Multi gene."""
    
    def __init__(self, input_dir, output_dir):
        """
        Initialize the alignment processor.
        
        Args:
            input_dir (str): Input directory path.
            output_dir (str): Output directory path.
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"入力ディレクトリ: {self.input_dir}")
        print(f"出力ディレクトリ: {self.output_dir}")
    
    def group_sequences_by_multi_id(self):
        """Group input sequences by Multi gene ID."""
        multi_groups = defaultdict(list)
        
        # Process all FASTA files in the input directory.
        for fasta_file in self.input_dir.glob("*.fa"):
            print(f"処理中: {fasta_file}")
            
            for record in SeqIO.parse(fasta_file, "fasta"):
                # Example: Multi_0000000001_GGGO01000001.1 -> Multi_0000000001
                header_parts = record.id.split('_')
                if len(header_parts) >= 2 and header_parts[0] == "Multi":
                    multi_id = f"{header_parts[0]}_{header_parts[1]}"
                    multi_groups[multi_id].append(record)
                    print(f"  {record.id} -> {multi_id}")
        
        print(f"見つかったMulti遺伝子グループ: {len(multi_groups)}")
        for multi_id, sequences in multi_groups.items():
            print(f"  {multi_id}: {len(sequences)} 配列")
        
        return multi_groups
    
    def run_mafft_alignment(self, sequences, multi_id):
        """Run MAFFT and return a validated multiple-sequence alignment."""
        if len(sequences) < 2:
            print(f"警告: {multi_id} に含まれる配列が1つのみです。アライメントをスキップします。")
            return sequences
        
        # Create temporary input and output files for MAFFT.
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as temp_input:
            temp_input_path = temp_input.name
            SeqIO.write(sequences, temp_input, "fasta")

        temp_output_path = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as temp_output:
                temp_output_path = temp_output.name
            
            cmd = [
                'mafft',
                '--auto',
                '--quiet',
                temp_input_path
            ]
            
            print(f"MAFFTアライメント実行中: {multi_id}")
            
            with open(temp_output_path, 'w') as output_file:
                result = subprocess.run(cmd, stdout=output_file, stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                error_message = (result.stderr or '').strip() or 'No error details provided.'
                raise RuntimeError(
                    f"MAFFT failed for {multi_id} with exit code "
                    f"{result.returncode}: {error_message}"
                )
            
            # Parse and validate the alignment before consensus generation.
            with open(temp_output_path, 'r') as alignment_file:
                aligned_sequences = list(SeqIO.parse(alignment_file, "fasta"))
            if not aligned_sequences:
                raise RuntimeError(
                    f"MAFFT produced no aligned sequences for {multi_id}."
                )
            if len(aligned_sequences) != len(sequences):
                raise RuntimeError(
                    f"MAFFT returned {len(aligned_sequences)} of "
                    f"{len(sequences)} sequences for {multi_id}."
                )

            alignment_lengths = {len(record.seq) for record in aligned_sequences}
            if len(alignment_lengths) != 1 or 0 in alignment_lengths:
                raise RuntimeError(
                    f"MAFFT produced an invalid alignment for {multi_id}: "
                    f"sequence lengths={sorted(alignment_lengths)}."
                )

            print(f"アライメント完了: {multi_id} ({len(aligned_sequences)} 配列)")
            
            return aligned_sequences
            
        finally:
            # Always clean up temporary files, including failure paths.
            if os.path.exists(temp_input_path):
                os.unlink(temp_input_path)
            if temp_output_path and os.path.exists(temp_output_path):
                os.unlink(temp_output_path)
    
    def generate_consensus(self, aligned_sequences, multi_id):
        """Generate a consensus sequence from a validated alignment."""
        if not aligned_sequences:
            return None
        
        if len(aligned_sequences) == 1:
            # A single sequence does not require alignment.
            consensus_seq = aligned_sequences[0].seq
            consensus_record = SeqRecord(consensus_seq, id=multi_id, description="consensus sequence")
            print(f"単一配列のコンセンサス: {multi_id}")
            return consensus_record
        
        alignment_lengths = {len(record.seq) for record in aligned_sequences}
        if len(alignment_lengths) != 1 or 0 in alignment_lengths:
            raise ValueError(
                f"Consensus input is not a valid alignment for {multi_id}: "
                f"sequence lengths={sorted(alignment_lengths)}."
            )

        # All records must have the same aligned length.
        alignment_length = len(aligned_sequences[0])
        
        # Select the most frequent non-gap base at each aligned position.
        consensus_sequence = []
        
        for pos in range(alignment_length):
            char_counts = defaultdict(int)
            
            for seq_record in aligned_sequences:
                if pos < len(seq_record.seq):
                    char = str(seq_record.seq[pos]).upper()
                    if char != '-':
                        char_counts[char] += 1
            
            if char_counts:
                consensus_char = max(char_counts, key=char_counts.get)
                consensus_sequence.append(consensus_char)
            else:
                consensus_sequence.append('N')
        
        # Remove alignment gaps from the final consensus.
        consensus_str = ''.join(consensus_sequence).replace('-', '')
        
        consensus_seq = Seq(consensus_str)
        consensus_record = SeqRecord(
            consensus_seq,
            id=multi_id,
            description=f"consensus sequence from {len(aligned_sequences)} sequences"
        )
        
        print(f"コンセンサス生成完了: {multi_id} (長さ: {len(consensus_str)})")
        return consensus_record
    
    def process_multi_gene_group(self, multi_id, sequences):
        """Align one Multi gene group and save its consensus sequence."""
        print(f"\n=== {multi_id} 処理開始 ===")
        
        aligned_sequences = self.run_mafft_alignment(sequences, multi_id)
        consensus_record = self.generate_consensus(aligned_sequences, multi_id)

        if consensus_record is None:
            raise RuntimeError(f"Consensus generation returned no result for {multi_id}.")

        output_filename = f"{multi_id}_aligned_consensus.fa"
        output_path = self.output_dir / output_filename
        SeqIO.write(consensus_record, output_path, "fasta")
        print(f"保存完了: {output_path}")
        
        print(f"=== {multi_id} 処理完了 ===")
        return output_path
    
    def run(self):
        """Process all groups and fail the run if any group fails."""
        print("Multi遺伝子アライメント処理を開始します...")
        
        multi_groups = self.group_sequences_by_multi_id()
        
        if not multi_groups:
            print("処理対象のMulti遺伝子グループが見つかりませんでした")
            return

        successful_count = 0
        failures = []

        for multi_id, sequences in multi_groups.items():
            try:
                self.process_multi_gene_group(multi_id, sequences)
                successful_count += 1
            except Exception as error:
                print(f"エラー: {multi_id} の処理中にエラーが発生しました: {error}")
                failures.append((multi_id, str(error)))
        
        print("\n=== 全体処理完了 ===")
        print(f"処理対象Multi遺伝子グループ: {len(multi_groups)}")
        print(f"成功したMulti遺伝子グループ: {successful_count}")
        print(f"失敗したMulti遺伝子グループ: {len(failures)}")
        
        output_files = list(self.output_dir.glob("*_aligned_consensus.fa"))
        print(f"生成されたファイル数: {len(output_files)}")

        if failures:
            failure_details = '; '.join(
                f"{multi_id}: {message}"
                for multi_id, message in failures
            )
            raise RuntimeError(
                f"Alignment failed for {len(failures)} of "
                f"{len(multi_groups)} groups: {failure_details}"
            )

        return {
            'total_groups': len(multi_groups),
            'successful_groups': successful_count,
            'failed_groups': 0,
        }


def main():
    """Run the alignment command-line interface."""
    parser = argparse.ArgumentParser(
        description='Align Multi-gene transcripts and generate consensus sequences.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python run_alignment_and_mask.py -i multi_genes -o multi_genes_aligned
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input directory path.'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory path.'
    )
    
    args = parser.parse_args()
    
    # Validate the input directory before creating output files.
    if not os.path.exists(args.input):
        print(f"エラー: 入力ディレクトリが存在しません: {args.input}")
        sys.exit(1)
    
    # Return a non-zero exit status when any group fails.
    processor = AlignmentProcessor(args.input, args.output)
    try:
        processor.run()
    except Exception as error:
        print(f"Alignment processing failed: {error}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
