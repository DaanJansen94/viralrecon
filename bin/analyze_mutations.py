#!/usr/bin/env python3

import os
import sys
import glob
import pandas as pd

def process_sample(variant_file, stats_file, sample_name):
    """Process a single sample to calculate mutation statistics."""
    
    # Get reads mapped from stats file
    with open(stats_file) as f:
        for line in f:
            if "reads mapped:" in line:
                reads = int(line.split()[3])
                break
    
    # Read variant file
    df = pd.read_csv(variant_file, sep='\t')
    
    # Filter for SNPs (where ALT is a single base)
    snp_mask = (~df['ALT'].str.startswith('-')) & (df['ALT'].str.len() == 1) & (df['ALT'].str.match('^[ATGC]$'))
    snps = df[snp_mask]
    
    # Calculate mutation counts
    all_muts = len(snps)
    major_muts = len(snps[snps['ALT_FREQ'] > 0.5])
    minor_muts = len(snps[snps['ALT_FREQ'] <= 0.5])
    
    # Calculate APOBEC3 mutations (G>A or C>T)
    apobec_mask = ((snps['REF'] == 'G') & (snps['ALT'] == 'A')) | ((snps['REF'] == 'C') & (snps['ALT'] == 'T'))
    major_apobec = len(snps[apobec_mask & (snps['ALT_FREQ'] > 0.5)])
    minor_apobec = len(snps[apobec_mask & (snps['ALT_FREQ'] <= 0.5)])
    
    # Calculate percentages
    major_apobec_perc = (major_apobec / major_muts * 100) if major_muts > 0 else 0.0
    minor_apobec_perc = (minor_apobec / minor_muts * 100) if minor_muts > 0 else 0.0
    
    # Calculate APM (APOBEC3 mutations Per Million MPXV reads)
    apm = (minor_apobec / reads * 1000000) if reads > 0 else 0.0
    
    return {
        'sample': sample_name,
        'reads': reads,
        'all_muts': all_muts,
        'major_muts': major_muts,
        'major_apobec': major_apobec,
        'major_apobec_perc': major_apobec_perc,
        'minor_muts': minor_muts,
        'minor_apobec': minor_apobec,
        'minor_apobec_perc': minor_apobec_perc,
        'apm': apm
    }

def main():
    if len(sys.argv) != 3:
        print("Usage: python analyze_mutations.py <variants_dir> <stats_dir>")
        print("Example: python analyze_mutations.py variants stats")
        sys.exit(1)
    
    variants_dir = sys.argv[1]
    stats_dir = sys.argv[2]
    output_file = "mutations_and_apobec3_summary.txt"
    
    # Check if input directories exist
    if not os.path.isdir(variants_dir):
        print(f"Error: Variants directory not found: {variants_dir}")
        sys.exit(1)
    if not os.path.isdir(stats_dir):
        print(f"Error: Stats directory not found: {stats_dir}")
        sys.exit(1)
    
    # Get variant files
    variant_files = glob.glob(os.path.join(variants_dir, "*.tsv"))
    print(f"Found variant files: {variant_files}")
    
    # Process all samples and collect results
    results = []
    for variant_file in variant_files:
        base_name = os.path.basename(variant_file)
        # Extract sample name without the .filtered.tsv suffix
        if '.filtered.tsv' in base_name:
            sample_name = base_name.replace('.filtered.tsv', '')
        else:
            sample_name = base_name.replace('.tsv', '')
            
        # Try to find the stats file with the correct naming pattern
        stats_file = os.path.join(stats_dir, f"{sample_name}.sorted.bam.stats")
        
        if not os.path.isfile(stats_file):
            print(f"Warning: Stats file not found for {base_name}: {stats_file}")
            continue
        
        print(f"Processing sample: {sample_name}")
        print(f"  Variant file: {variant_file}")
        print(f"  Stats file: {stats_file}")
        
        result = process_sample(variant_file, stats_file, sample_name)
        results.append(result)
    
    if not results:
        print("Error: No samples could be processed")
        sys.exit(1)
    
    # Calculate totals
    total_reads = sum(r['reads'] for r in results)
    total_muts = sum(r['all_muts'] for r in results)
    total_major = sum(r['major_muts'] for r in results)
    total_major_apobec = sum(r['major_apobec'] for r in results)
    total_minor = sum(r['minor_muts'] for r in results)
    total_minor_apobec = sum(r['minor_apobec'] for r in results)
    
    total_major_perc = (total_major_apobec / total_major * 100) if total_major > 0 else 0.0
    total_minor_perc = (total_minor_apobec / total_minor * 100) if total_minor > 0 else 0.0
    total_apm = (total_minor_apobec / total_reads * 1000000) if total_reads > 0 else 0.0
    
    # Write output
    with open(output_file, 'w') as f:
        f.write("Analysis of mutations (SNPs only, excluding insertions and deletions)\n")
        f.write("================================\n\n")
        f.write("NOTE: This analysis focuses on single nucleotide polymorphisms (SNPs) only.\n")
        f.write("Results do not include insertions, deletions, or other complex variants.\n")
        f.write("Major mutations: frequency > 50%\n")
        f.write("Minor mutations: frequency <= 50%\n")
        f.write("APM: APOBEC3 mutations Per Million MPXV reads\n\n")
        f.write("Summary Table - Mutations by Sample\n")
        f.write("------------------------------------------\n\n")
        
        # Write header
        f.write(f"{'Sample':<30} {'Reads':<10} {'All Muts':<10} {'Major Muts':<10} {'Major APOB':<10} {'Major %':<10} {'Minor Muts':<10} {'Minor APOB':<10} {'Minor %':<10} {'APM':<10}\n")
        f.write('-' * 130 + '\n')
        
        # Write each sample
        for r in results:
            f.write(f"{r['sample']:<30} {r['reads']:<10} {r['all_muts']:<10} {r['major_muts']:<10} {r['major_apobec']:<10} {r['major_apobec_perc']:<10.2f} {r['minor_muts']:<10} {r['minor_apobec']:<10} {r['minor_apobec_perc']:<10.2f} {r['apm']:<10.2f}\n")
        
        # Write totals
        f.write('-' * 130 + '\n')
        f.write(f"{'TOTAL':<30} {total_reads:<10} {total_muts:<10} {total_major:<10} {total_major_apobec:<10} {total_major_perc:<10.2f} {total_minor:<10} {total_minor_apobec:<10} {total_minor_perc:<10.2f} {total_apm:<10.2f}\n")
    
    print(f"Results written to {output_file}")
    
if __name__ == '__main__':
    main() 