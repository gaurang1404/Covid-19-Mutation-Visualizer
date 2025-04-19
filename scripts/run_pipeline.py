#!/usr/bin/env python
# run_pipeline.py

import os
import sys
import argparse
import time
from enhanced_alignment import align_sequences
from mutation_analysis import analyze_mutations

def main():
    parser = argparse.ArgumentParser(description='COVID-19 Mutation Tracker and Analyzer')
    parser.add_argument('--ref', type=str, default='data/ref.fasta', help='Path to reference FASTA file')
    parser.add_argument('--samples', type=str, nargs='+', help='Paths to sample FASTA files')
    parser.add_argument('--output', type=str, default='results', help='Output directory')
    
    args = parser.parse_args()
    
    # If no sample files specified, look for them in the data directory
    if not args.samples:
        data_dir = os.path.dirname(args.ref)
        samples = [os.path.join(data_dir, f) for f in os.listdir(data_dir) 
                  if f.endswith('.fasta') and not f == os.path.basename(args.ref)]
        
        if not samples:
            print("Error: No sample FASTA files found in data directory.")
            sys.exit(1)
        
        args.samples = samples
    
    # Create output directories
    os.makedirs(args.output, exist_ok=True)
    visualization_dir = os.path.join(os.path.dirname(args.output), 'visualizations')
    os.makedirs(visualization_dir, exist_ok=True)
    
    print("COVID-19 Mutation Tracker and Analyzer")
    print("====================================")
    print(f"Reference: {args.ref}")
    print(f"Samples: {', '.join(args.samples)}")
    print(f"Output directory: {args.output}")
    print("\nStarting analysis pipeline...")
    
    start_time = time.time()
    
    # Step 1: Align sequences and detect mutations
    print("\n1. Performing sequence alignment and mutation detection...")
    mutations_df = align_sequences(args.ref, args.samples, args.output)
    
    # Step 2: Analyze mutations if any were found
    if not mutations_df.empty:
        print("\n2. Analyzing mutation patterns...")
        mutations_file = os.path.join(args.output, "all_mutations.csv")
        analyze_mutations(mutations_file, visualization_dir)
    else:
        print("\nNo mutations were detected. Analysis skipped.")
    
    elapsed_time = time.time() - start_time
    print(f"\nAnalysis complete in {elapsed_time:.2f} seconds")
    print(f"Results saved to {args.output} and {visualization_dir}")

if __name__ == "__main__":
    main()