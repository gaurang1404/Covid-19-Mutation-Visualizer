#!/usr/bin/env python
# alignment.py

import os
from Bio import SeqIO, Align, Seq

import pandas as pd
from collections import defaultdict

def load_fasta(fasta_file):
    """Load sequences from a FASTA file"""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def perform_pairwise_alignment(ref_seq, query_seq):
    """Perform local alignment between reference and query sequences"""
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'  # Use local alignment to handle sequences of different lengths
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    
    # Perform alignment
    alignments = aligner.align(ref_seq, query_seq)
    # Return the best alignment
    return alignments[0]

def find_mutations_from_alignment(alignment, ref_seq, query_seq, ref_start=0):
    """Extract mutations from an alignment object"""
    mutations = []
    
    # Convert alignment to string representation
    alignment_str = str(alignment)
    alignment_lines = alignment_str.split('\n')
    
    if len(alignment_lines) < 3:
        print("Warning: Unexpected alignment format")
        return mutations
    
    # Get the aligned sequences
    ref_aligned = alignment_lines[0]
    query_aligned = alignment_lines[2]
    
    # Determine the start position in the reference
    ref_pos = ref_start
    query_pos = 0
    
    # Parse the alignment to find mutations
    for i in range(len(ref_aligned)):
        ref_char = ref_aligned[i]
        query_char = query_aligned[i]
        
        # Skip alignment spacers
        if ref_char == ' ' or query_char == ' ':
            continue
        
        # Track positions
        if ref_char != '-':
            ref_pos += 1
        if query_char != '-':
            query_pos += 1
        
        # Check for mutations
        if ref_char != query_char:
            if ref_char == '-':
                # Insertion
                # Find the complete inserted sequence
                ins_start = i
                ins_end = i
                while ins_end < len(ref_aligned) and ref_aligned[ins_end] == '-':
                    ins_end += 1
                inserted_seq = query_aligned[ins_start:ins_end]
                
                mutations.append({
                    'Position': ref_pos,
                    'Ref': '-',
                    'Alt': inserted_seq,
                    'Type': 'INS'
                })
            elif query_char == '-':
                # Deletion
                mutations.append({
                    'Position': ref_pos,
                    'Ref': ref_char,
                    'Alt': '-',
                    'Type': 'DEL'
                })
            else:
                # SNP (Single Nucleotide Polymorphism)
                mutations.append({
                    'Position': ref_pos,
                    'Ref': ref_char,
                    'Alt': query_char,
                    'Type': 'SNP'
                })
    
    return mutations

def process_large_sequences(ref_seq, query_seq, window_size=10000, overlap=1000):
    """Handle large sequences by processing them in windows"""
    all_mutations = []
    ref_length = len(ref_seq)
    
    # Process the sequence in overlapping windows
    for i in range(0, ref_length, window_size - overlap):
        window_end = min(i + window_size, ref_length)
        ref_window = ref_seq[i:window_end]
        
        # Try to find the best matching region in the query sequence
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -5
        aligner.extend_gap_score = -0.5
        
        # Align the reference window to the query
        alignments = aligner.align(ref_window, query_seq)
        best_alignment = alignments[0]
        
        # Extract mutations from this window's alignment
        window_mutations = find_mutations_from_alignment(best_alignment, ref_window, query_seq, ref_start=i)
        all_mutations.extend(window_mutations)
    
    # Remove potential duplicate mutations from overlapping regions
    unique_mutations = {}
    for mut in all_mutations:
        key = f"{mut['Position']}:{mut['Ref']}>{mut['Alt']}"
        unique_mutations[key] = mut
    
    return list(unique_mutations.values())

def align_sequences(reference_file, sample_files, output_dir):
    """Align sample sequences to the reference using Biopython"""
    
    # Load reference sequence
    ref_sequences = load_fasta(reference_file)
    reference_id = list(ref_sequences.keys())[0]
    reference_seq = ref_sequences[reference_id]
    
    print(f"Reference sequence: {reference_id}, Length: {len(reference_seq)}")
    
    mutation_data = []
    
    # Process each sample file
    for sample_file in sample_files:
        country = os.path.basename(sample_file).split('.')[0]
        print(f"\nProcessing sample from {country}")
        
        sample_sequences = load_fasta(sample_file)
        
        for seq_id, sequence in sample_sequences.items():
            print(f"  Sequence: {seq_id}, Length: {len(sequence)}")
            
            # Choose alignment strategy based on sequence lengths
            if len(reference_seq) > 50000 or len(sequence) > 50000:
                print("  Using windowed alignment approach for large sequences")
                mutations = process_large_sequences(reference_seq, sequence)
            else:
                # For smaller sequences, do a direct alignment
                alignment = perform_pairwise_alignment(reference_seq, sequence)
                mutations = find_mutations_from_alignment(alignment, reference_seq, sequence)
            
            print(f"  Found {len(mutations)} mutations")
            
            # Add country and sample ID to each mutation
            for mut in mutations:
                mut['Country'] = country
                mut['Sample_ID'] = seq_id
                mutation_data.append(mut)
    
    # Create mutations DataFrame and save to CSV
    if mutation_data:
        mutations_df = pd.DataFrame(mutation_data)
        
        # Save complete mutation data
        output_file = os.path.join(output_dir, "all_mutations.csv")
        mutations_df.to_csv(output_file, index=False)
        print(f"\nMutation data saved to {output_file}")
        
        # Summary by country
        country_summary = mutations_df.groupby(['Country', 'Type']).size().unstack(fill_value=0)
        summary_file = os.path.join(output_dir, "mutation_summary_by_country.csv")
        country_summary.to_csv(summary_file)
        print(f"Country summary saved to {summary_file}")
        
        return mutations_df
    else:
        print("No mutations detected")
        return pd.DataFrame()

if __name__ == "__main__":
    # Paths to reference and sample files
    reference_file = "data/ref.fasta"
    sample_files = [f"data/{country}.fasta" for country in ["sample1", "sample2", "sample3", "sample4", "sample5"]]
    output_dir = "results"
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Run alignment and mutation detection
    align_sequences(reference_file, sample_files, output_dir)