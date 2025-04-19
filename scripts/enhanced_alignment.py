#!/usr/bin/env python
# enhanced_alignment.py

import os

from Bio import SeqIO, Seq, Align
from Bio.Align import substitution_matrices

import pandas as pd
import numpy as np
from collections import defaultdict
import multiprocessing as mp
from functools import partial

def load_fasta(fasta_file):
    """Load sequences from a FASTA file"""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def perform_local_alignment(ref_seq, query_seq):
    """Perform local alignment between reference and query sequences"""
    # Create aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    # Perform alignment
    alignments = aligner.align(ref_seq, query_seq)
    
    # Return the best alignment
    if len(alignments) > 0:
        return alignments[0]
    else:
        return None

def extract_mutations_from_alignment(alignment, country, sample_id):
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
    
    # Get alignment coordinates - in the new API we can use these properties
    start = alignment.coordinates[0][0]  # Start position in reference sequence
    
    # Process the alignment and identify mutations
    ref_pos = start
    
    for i in range(len(ref_aligned)):
        ref_char = ref_aligned[i]
        query_char = query_aligned[i]
        
        # Skip alignment spacers
        if ref_char == ' ' or query_char == ' ':
            continue
        
        # Handle insertions (gaps in reference)
        if ref_char == '-':
            # Skip counting reference position
            pass
        # Handle deletions (gaps in query)
        elif query_char == '-':
            ref_pos += 1
        # Handle matches and mismatches
        else:
            if ref_char != query_char:
                mutations.append({
                    'Country': country,
                    'Sample_ID': sample_id,
                    'Position': ref_pos + 1,  # 1-based position
                    'Ref': ref_char,
                    'Alt': query_char,
                    'Type': 'SNP'
                })
            ref_pos += 1
    
    return mutations

def process_sequence_chunk(chunk_data):
    """Process a chunk of the reference sequence against a query sequence"""
    ref_chunk, ref_start, query_seq, country, sample_id = chunk_data
    
    # Perform local alignment on this chunk
    alignment = perform_local_alignment(ref_chunk, query_seq)
    
    if alignment:
        # Extract mutations, adjusting positions to reflect global coordinates
        local_mutations = extract_mutations_from_alignment(alignment, country, sample_id)
        
        # Adjust positions by adding the chunk start position
        for mut in local_mutations:
            mut['Position'] += ref_start
        
        return local_mutations
    else:
        return []

def align_sequences_chunked(reference_file, sample_files, output_dir, chunk_size=2000, overlap=200):
    """Align sequences using a chunking approach to handle different lengths efficiently"""
    
    # Load reference sequence
    ref_sequences = load_fasta(reference_file)
    reference_id = list(ref_sequences.keys())[0]
    reference_seq = ref_sequences[reference_id]
    
    print(f"Reference sequence: {reference_id}, Length: {len(reference_seq)}")
    
    all_mutations = []
    
    # Process each sample file
    for sample_file in sample_files:
        country = os.path.basename(sample_file).split('.')[0]
        print(f"\nProcessing sample from {country}")
        
        sample_sequences = load_fasta(sample_file)
        
        for seq_id, query_seq in sample_sequences.items():
            print(f"  Sequence: {seq_id}, Length: {len(query_seq)}")
            
            # Break reference into overlapping chunks
            chunks = []
            
            for i in range(0, len(reference_seq), chunk_size - overlap):
                chunk_end = min(i + chunk_size, len(reference_seq))
                ref_chunk = reference_seq[i:chunk_end]
                chunks.append((ref_chunk, i, query_seq, country, seq_id))
            
            print(f"  Processing in {len(chunks)} chunks")
            
            # Process chunks in parallel
            with mp.Pool(processes=mp.cpu_count()) as pool:
                chunk_results = pool.map(process_sequence_chunk, chunks)
            
            # Combine results
            sequence_mutations = []
            for result in chunk_results:
                sequence_mutations.extend(result)
            
            # Remove duplicate mutations (might occur in overlapping regions)
            unique_mutations = {}
            for mut in sequence_mutations:
                key = f"{mut['Position']}:{mut['Ref']}>{mut['Alt']}"
                unique_mutations[key] = mut
            
            sequence_mutations = list(unique_mutations.values())
            all_mutations.extend(sequence_mutations)
            
            print(f"  Found {len(sequence_mutations)} unique mutations")
    
    # Create mutations DataFrame and save to CSV
    if all_mutations:
        mutations_df = pd.DataFrame(all_mutations)
        
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

#!/usr/bin/env python
# enhanced_alignment.py

import os
from Bio import SeqIO, Seq, Align
from Bio.Align import substitution_matrices
import pandas as pd
import numpy as np
from collections import defaultdict
import multiprocessing as mp
from functools import partial

def load_fasta(fasta_file):
    """Load sequences from a FASTA file"""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def perform_local_alignment(ref_seq, query_seq):
    """Perform local alignment between reference and query sequences using Bio.Align"""
    # Create aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    
    # Perform alignment
    alignments = aligner.align(ref_seq, query_seq)
    
    # Return the best alignment
    if len(alignments) > 0:
        return alignments[0]
    else:
        return None

def extract_mutations_from_alignment(alignment, country, sample_id):
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
    
    # Get alignment coordinates
    start = alignment.coordinates[0][0]  # Start position in reference sequence
    
    # Process the alignment and identify mutations
    ref_pos = start
    
    for i in range(len(ref_aligned)):
        ref_char = ref_aligned[i]
        query_char = query_aligned[i]
        
        # Skip alignment spacers
        if ref_char == ' ' or query_char == ' ':
            continue
        
        # Handle insertions (gaps in reference)
        if ref_char == '-':
            # Skip counting reference position
            pass
        # Handle deletions (gaps in query)
        elif query_char == '-':
            ref_pos += 1
        # Handle matches and mismatches
        else:
            if ref_char != query_char:
                mutations.append({
                    'Country': country,
                    'Sample_ID': sample_id,
                    'Position': ref_pos + 1,  # 1-based position
                    'Ref': ref_char,
                    'Alt': query_char,
                    'Type': 'SNP'
                })
            ref_pos += 1
    
    return mutations

def process_sequence_chunk(chunk_data):
    """Process a chunk of the reference sequence against a query sequence"""
    ref_chunk, ref_start, query_seq, country, sample_id = chunk_data
    
    # Perform local alignment on this chunk
    alignment = perform_local_alignment(ref_chunk, query_seq)
    
    if alignment:
        # Extract mutations, adjusting positions to reflect global coordinates
        local_mutations = extract_mutations_from_alignment(alignment, country, sample_id)
        
        # Adjust positions by adding the chunk start position
        for mut in local_mutations:
            mut['Position'] += ref_start
        
        return local_mutations
    else:
        return []

def align_sequences_chunked(reference_file, sample_files, output_dir, chunk_size=2000, overlap=200):
    """Align sequences using a chunking approach to handle different lengths efficiently"""
    
    # Load reference sequence
    ref_sequences = load_fasta(reference_file)
    reference_id = list(ref_sequences.keys())[0]
    reference_seq = ref_sequences[reference_id]
    
    print(f"Reference sequence: {reference_id}, Length: {len(reference_seq)}")
    
    all_mutations = []
    
    # Process each sample file
    for sample_file in sample_files:
        country = os.path.basename(sample_file).split('.')[0]
        print(f"\nProcessing sample from {country}")
        
        sample_sequences = load_fasta(sample_file)
        
        for seq_id, query_seq in sample_sequences.items():
            print(f"  Sequence: {seq_id}, Length: {len(query_seq)}")
            
            # Break reference into overlapping chunks
            chunks = []
            
            for i in range(0, len(reference_seq), chunk_size - overlap):
                chunk_end = min(i + chunk_size, len(reference_seq))
                ref_chunk = reference_seq[i:chunk_end]
                chunks.append((ref_chunk, i, query_seq, country, seq_id))
            
            print(f"  Processing in {len(chunks)} chunks")
            
            # Process chunks in parallel
            with mp.Pool(processes=mp.cpu_count()) as pool:
                chunk_results = pool.map(process_sequence_chunk, chunks)
            
            # Combine results
            sequence_mutations = []
            for result in chunk_results:
                sequence_mutations.extend(result)
            
            # Remove duplicate mutations (might occur in overlapping regions)
            unique_mutations = {}
            for mut in sequence_mutations:
                key = f"{mut['Position']}:{mut['Ref']}>{mut['Alt']}"
                unique_mutations[key] = mut
            
            sequence_mutations = list(unique_mutations.values())
            all_mutations.extend(sequence_mutations)
            
            print(f"  Found {len(sequence_mutations)} unique mutations")
    
    # Create mutations DataFrame and save to CSV
    if all_mutations:
        mutations_df = pd.DataFrame(all_mutations)
        
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

def detect_indels(reference_file, sample_files, output_dir, window_size=100):
    """Detect insertions and deletions by comparing sequence lengths in sliding windows"""
    
    ref_sequences = load_fasta(reference_file)
    reference_id = list(ref_sequences.keys())[0]
    reference_seq = ref_sequences[reference_id]
    
    indel_mutations = []
    
    for sample_file in sample_files:
        country = os.path.basename(sample_file).split('.')[0]
        print(f"\nDetecting indels in sample from {country}")
        
        sample_sequences = load_fasta(sample_file)
        
        for seq_id, query_seq in sample_sequences.items():
            print(f"  Sequence: {seq_id}")
            
            # Use dynamic programming to find optimal local alignments
            # that might contain indels
            for i in range(0, len(reference_seq) - window_size, window_size // 2):
                ref_window = reference_seq[i:i+window_size]
                
                # Find best matching region in query sequence
                best_score = -1
                best_j = -1
                
                # Simple sliding window search in query
                for j in range(0, max(len(query_seq) - window_size, 1), 10):  # Step by 10 for efficiency
                    query_window = query_seq[j:j+window_size]
                    
                    # Count matching bases (simple scoring)
                    score = sum(1 for r, q in zip(ref_window, query_window) if r == q)
                    
                    if score > best_score:
                        best_score = score
                        best_j = j
                
                if best_j >= 0:
                    # Find exact alignment with indels
                    aligner = Align.PairwiseAligner()
                    aligner.mode = 'local'
                    aligner.match_score = 2
                    aligner.mismatch_score = -1
                    aligner.open_gap_score = -7
                    aligner.extend_gap_score = -0.5
                    
                    alignments = aligner.align(ref_window, query_seq[best_j:best_j+window_size+50])
                    
                    if len(alignments) > 0:
                        alignment_obj = alignments[0]
                        alignment_str = str(alignment_obj)
                        alignment_lines = alignment_str.split('\n')
                        
                        if len(alignment_lines) >= 3:
                            ref_aligned = alignment_lines[0]
                            query_aligned = alignment_lines[2]
                            
                            # Scan alignment for indels
                            ref_pos = i
                            indel_start = None
                            indel_type = None
                            indel_ref = ""
                            indel_alt = ""
                            
                            for k in range(len(ref_aligned)):
                                if k < len(query_aligned):  # Make sure we don't go out of bounds
                                    if ref_aligned[k] == '-' and indel_type != 'INS':
                                        # Start of insertion or continue insertion
                                        if indel_type is None:
                                            indel_start = ref_pos
                                            indel_type = 'INS'
                                            indel_ref = '-'
                                            indel_alt = query_aligned[k]
                                        else:
                                            # Close previous indel
                                            indel_mutations.append({
                                                'Country': country,
                                                'Sample_ID': seq_id,
                                                'Position': indel_start + 1,  # 1-based position
                                                'Ref': indel_ref,
                                                'Alt': indel_alt,
                                                'Type': indel_type
                                            })
                                            # Start new indel
                                            indel_start = ref_pos
                                            indel_type = 'INS'
                                            indel_ref = '-'
                                            indel_alt = query_aligned[k]
                                    elif query_aligned[k] == '-' and indel_type != 'DEL':
                                        # Start of deletion or continue deletion
                                        if indel_type is None:
                                            indel_start = ref_pos
                                            indel_type = 'DEL'
                                            indel_ref = ref_aligned[k]
                                            indel_alt = '-'
                                        else:
                                            # Close previous indel
                                            indel_mutations.append({
                                                'Country': country,
                                                'Sample_ID': seq_id,
                                                'Position': indel_start + 1,  # 1-based position
                                                'Ref': indel_ref,
                                                'Alt': indel_alt,
                                                'Type': indel_type
                                            })
                                            # Start new indel
                                            indel_start = ref_pos
                                            indel_type = 'DEL'
                                            indel_ref = ref_aligned[k]
                                            indel_alt = '-'
                                    elif ref_aligned[k] == '-' and indel_type == 'INS':
                                        # Continue insertion
                                        indel_alt += query_aligned[k]
                                    elif query_aligned[k] == '-' and indel_type == 'DEL':
                                        # Continue deletion
                                        indel_ref += ref_aligned[k]
                                    else:
                                        # End of indel if there was one
                                        if indel_type is not None:
                                            indel_mutations.append({
                                                'Country': country,
                                                'Sample_ID': seq_id,
                                                'Position': indel_start + 1,  # 1-based position
                                                'Ref': indel_ref,
                                                'Alt': indel_alt,
                                                'Type': indel_type
                                            })
                                            indel_type = None
                                        
                                        # Increment reference position for matches/mismatches
                                        if ref_aligned[k] != '-':
                                            ref_pos += 1
    
    # Save indel mutations
    if indel_mutations:
        indels_df = pd.DataFrame(indel_mutations)
        output_file = os.path.join(output_dir, "indel_mutations.csv")
        indels_df.to_csv(output_file, index=False)
        print(f"\nIndel mutation data saved to {output_file}")
        
        return indels_df
    else:
        print("No indels detected")
        return pd.DataFrame()

def align_sequences(reference_file, sample_files, output_dir):
    """Main function to align sequences and detect all mutation types"""
    
    # First perform chunked alignment to detect SNPs
    snp_df = align_sequences_chunked(reference_file, sample_files, output_dir)
    
    # Then detect indels
    indel_df = detect_indels(reference_file, sample_files, output_dir)
    
    # Combine SNPs and indels
    if not snp_df.empty and not indel_df.empty:
        combined_df = pd.concat([snp_df, indel_df])
        
        # Save combined mutations
        output_file = os.path.join(output_dir, "all_mutations.csv")
        combined_df.to_csv(output_file, index=False)
        print(f"\nCombined mutation data saved to {output_file}")
        
        return combined_df
    elif not snp_df.empty:
        return snp_df
    elif not indel_df.empty:
        return indel_df
    else:
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