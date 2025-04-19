#!/usr/bin/env python
# mutation_analysis.py

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def analyze_mutations(mutations_file, output_dir):
    """Analyze mutations and create summary statistics"""
    
    # Load mutation data
    df = pd.read_csv(mutations_file)
    print(f"Loaded {len(df)} mutations")
    
    # 1. Mutation type distribution by country
    plt.figure(figsize=(12, 6))
    country_type_counts = df.groupby(['Country', 'Type']).size().unstack(fill_value=0)
    country_type_counts.plot(kind='bar', stacked=True)
    plt.title('Mutation Types by Country')
    plt.xlabel('Country')
    plt.ylabel('Number of Mutations')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'mutation_types_by_country.png'))
    
    # 2. Mutation hotspots - find regions with high mutation rates
    plt.figure(figsize=(14, 6))
    position_counts = df.groupby('Position').size()
    position_counts.plot(kind='line')
    plt.title('Mutation Hotspots Across the Genome')
    plt.xlabel('Genome Position')
    plt.ylabel('Mutation Count')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'mutation_hotspots.png'))
    
    # 3. SNP substitution patterns
    if 'SNP' in df['Type'].values:
        snps = df[df['Type'] == 'SNP']
        substitutions = snps.groupby(['Ref', 'Alt']).size().unstack(fill_value=0)
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(substitutions, annot=True, cmap='YlGnBu', fmt='d')
        plt.title('SNP Substitution Patterns')
        plt.savefig(os.path.join(output_dir, 'snp_substitution_patterns.png'))
    
    # 4. Shared mutations across countries
    mutation_keys = df.apply(lambda x: f"{x['Position']}:{x['Ref']}>{x['Alt']}", axis=1)
    df['Mutation'] = mutation_keys
    
    shared_mutations = defaultdict(list)
    for mutation in df['Mutation'].unique():
        countries = df[df['Mutation'] == mutation]['Country'].unique()
        shared_mutations[tuple(sorted(countries))].append(mutation)
    
    # Create a summary of shared mutations
    shared_summary = []
    for countries, mutations in shared_mutations.items():
        if len(countries) > 1:  # Only include mutations shared across multiple countries
            shared_summary.append({
                'Countries': ', '.join(countries),
                'Num_Countries': len(countries),
                'Shared_Mutations': len(mutations),
                'Mutations': '; '.join(mutations[:5]) + ('...' if len(mutations) > 5 else '')
            })
    
    if shared_summary:
        shared_df = pd.DataFrame(shared_summary)
        shared_df = shared_df.sort_values('Num_Countries', ascending=False)
        shared_df.to_csv(os.path.join(output_dir, 'shared_mutations.csv'), index=False)
        
        # Visualize mutations shared across countries
        plt.figure(figsize=(12, 6))
        top_shared = shared_df.head(10)  # Top 10 shared mutation patterns
        sns.barplot(x='Num_Countries', y='Countries', data=top_shared, 
                   hue='Shared_Mutations', palette='viridis')
        plt.title('Top Mutation Patterns Shared Across Countries')
        plt.xlabel('Number of Countries')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'shared_mutation_patterns.png'))
    
    # 5. Create a comprehensive mutation report
    with open(os.path.join(output_dir, 'mutation_report.txt'), 'w') as f:
        f.write("COVID-19 Mutation Analysis Report\n")
        f.write("===============================\n\n")
        
        f.write(f"Total mutations detected: {len(df)}\n")
        f.write(f"Mutation types: {df['Type'].value_counts().to_dict()}\n\n")
        
        f.write("Mutations by country:\n")
        country_counts = df.groupby('Country').size()
        for country, count in country_counts.items():
            f.write(f"  {country}: {count} mutations\n")
        
        if 'SNP' in df['Type'].values:
            f.write("\nSNP analysis:\n")
            f.write(f"  Total SNPs: {len(snps)}\n")
            
            # Most frequent SNP locations
            top_snp_positions = snps['Position'].value_counts().head(5)
            f.write("  Top 5 SNP positions:\n")
            for pos, count in top_snp_positions.items():
                f.write(f"    Position {pos}: {count} occurrences\n")
        
        if shared_summary:
            f.write("\nShared mutations analysis:\n")
            f.write(f"  {len(shared_summary)} unique mutation patterns shared across multiple countries\n")
            f.write("  Top 3 widely shared mutation patterns:\n")
            for _, row in shared_df.head(3).iterrows():
                f.write(f"    Shared across {row['Num_Countries']} countries ({row['Countries']}): {row['Shared_Mutations']} mutations\n")
    
    print("Analysis complete. Results saved to the visualizations directory.")

if __name__ == "__main__":
    mutations_file = "results/all_mutations.csv"
    output_dir = "visualizations"
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Run analysis
    analyze_mutations(mutations_file, output_dir)