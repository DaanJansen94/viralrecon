#!/usr/bin/env python

import pysam
import pandas as pd
import os
import glob
import sys
import shutil

def main():
    # Get command line arguments
    if len(sys.argv) != 3:
        print("Usage: apobec3_analysis.py <vcf_file> <reference_file>")
        sys.exit(1)
        
    vcf_file = sys.argv[1]
    ref_file = sys.argv[2]

    # Initialize storage for all results
    all_apobec_mutations = []

    print(f"Processing {vcf_file}...")
    
    # Get sample name from VCF file (remove .gz extension if present)
    sample_name = os.path.splitext(vcf_file)[0]
    if sample_name.endswith('.vcf'):
        sample_name = os.path.splitext(sample_name)[0]
    
    # Open reference genome and VCF
    try:
        ref = pysam.FastaFile(ref_file)
        vcf = pysam.VariantFile(vcf_file)
    except Exception as e:
        print(f"Error opening files for {vcf_file}: {e}")
        sys.exit(1)

    # Initialize storage for this sample's mutations
    sample_mutations = []
    tc_tt_mutations = []
    ga_aa_mutations = []

    # Iterate through variants in the VCF file
    for variant in vcf:
        chrom = variant.chrom
        pos = variant.pos  # 1-based position
        ref_allele = variant.ref
        
        # Skip if no alternate alleles or if it's not a SNP
        if not variant.alts or len(ref_allele) != 1 or len(variant.alts[0]) != 1:
            continue
            
        alt_allele = variant.alts[0]
        
        try:
            # Get information common to all mutations
            mutation_info = {
                'chrom': chrom,
                'pos': pos,
                'ref': ref_allele,
                'alt': alt_allele,
                'depth': variant.info.get('DP', 0) if 'DP' in variant.info else 0,
                'quality': variant.qual,
                'sample': sample_name
            }
            
            # Check for TC>TT (C>T where C is preceded by T)
            if pos > 1:  # Ensure we're not at the start of chromosome
                preceding_base = ref.fetch(chrom, pos-2, pos-1).upper()
                dinuc_ref = preceding_base + ref_allele
                dinuc_alt = preceding_base + alt_allele
                
                if dinuc_ref == "TC" and dinuc_alt == "TT":
                    mutation_info_copy = mutation_info.copy()
                    mutation_info_copy['context'] = dinuc_ref + ">" + dinuc_alt
                    mutation_info_copy['type'] = "TC>TT"
                    tc_tt_mutations.append(mutation_info_copy)
                    sample_mutations.append(mutation_info_copy)
                    all_apobec_mutations.append(mutation_info_copy)
            
            # Check for GA>AA (G>A where G is followed by A)
            following_base = ref.fetch(chrom, pos, pos+1).upper()
            dinuc_ref_rev = ref_allele + following_base
            dinuc_alt_rev = alt_allele + following_base
            
            if dinuc_ref_rev == "GA" and dinuc_alt_rev == "AA":
                mutation_info_copy = mutation_info.copy()
                mutation_info_copy['context'] = dinuc_ref_rev + ">" + dinuc_alt_rev
                mutation_info_copy['type'] = "GA>AA"
                ga_aa_mutations.append(mutation_info_copy)
                sample_mutations.append(mutation_info_copy)
                all_apobec_mutations.append(mutation_info_copy)
                
        except Exception as e:
            print(f"Error processing variant at {chrom}:{pos}: {e}")

    # Create DataFrame for this sample and save to CSV
    sample_df = pd.DataFrame(sample_mutations)
    if not sample_df.empty:
        sample_output_file = f"{sample_name}_apobec_mutations.csv"
        
        # Calculate sample summary statistics
        # Ensure total mutations equals the sum of TC>TT and GA>AA
        tc_tt_count = len(tc_tt_mutations)
        ga_aa_count = len(ga_aa_mutations)
        total_mutations = tc_tt_count + ga_aa_count  # Should equal len(sample_mutations)
        
        # Verify the counts match
        if len(sample_mutations) != total_mutations:
            print(f"Warning: Count mismatch - Sample mutations: {len(sample_mutations)}, TC>TT + GA>AA: {total_mutations}")
            print(f"Adjusting total count to match sum of subtypes")
        
        # Save sample mutations to CSV
        sample_df.to_csv(sample_output_file, index=False)
        
        # Append summary statistics to the sample CSV file
        with open(sample_output_file, 'a') as f:
            f.write('\n\nSummary Statistics\n')
            f.write(f'Total APOBEC3-related mutations,{total_mutations}\n')
            f.write(f'TC>TT mutations,{tc_tt_count}\n')
            f.write(f'GA>AA mutations,{ga_aa_count}\n')
        
        print(f"Sample mutations saved to {sample_output_file}")
    
    # Print summary for this sample
    print(f"\nSummary for {vcf_file}:")
    print(f"Total TC>TT mutations found: {tc_tt_count}")
    print(f"Total GA>AA mutations found: {ga_aa_count}")
    print(f"Total APOBEC3-related mutations: {total_mutations}")

    # Create a complete summary file
    # Instead of just finding files in the current directory,
    # Save current sample data to a dictionary
    sample_data = {}
    
    # First check for files in the current directory (which includes the one we just processed)
    current_files = glob.glob("*_apobec_mutations.csv")
    print(f"Found {len(current_files)} sample file(s) in current directory")
    
    # Extract data from current sample
    for file in current_files:
        try:
            file_sample = os.path.basename(file).replace("_apobec_mutations.csv", "")
            df = pd.read_csv(file)
            
            # Count mutations by type
            tc_tt = len(df[df['type'] == 'TC>TT']) if 'type' in df.columns else 0
            ga_aa = len(df[df['type'] == 'GA>AA']) if 'type' in df.columns else 0
            
            # Total should be the sum of the two types
            total_count = tc_tt + ga_aa
            
            sample_data[file_sample] = {
                'total': total_count,
                'tc_tt': tc_tt,
                'ga_aa': ga_aa
            }
            print(f"Added sample {file_sample} with {total_count} mutations ({tc_tt} TC>TT, {ga_aa} GA>AA)")
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    # Store basic summary data for current sample in a file that can be combined later
    current_sample_summary = {
        'sample': sample_name,
        'total': total_mutations,
        'tc_tt': tc_tt_count,
        'ga_aa': ga_aa_count
    }
    
    # Write this data to a special file that can be used for combining later
    with open(f"{sample_name}_summary_data.txt", 'w') as f:
        f.write(f"Sample: {sample_name}\n")
        f.write(f"Total: {total_mutations}\n")
        f.write(f"TC>TT: {tc_tt_count}\n")
        f.write(f"GA>AA: {ga_aa_count}\n")
    
    # Write summary file
    with open("apobec_summary.txt", 'w') as f:
        f.write("APOBEC3 Mutation Analysis Summary\n")
        f.write("================================\n\n")
        
        # Create the summary table at the top
        f.write("Summary Table - APOBEC3 Mutations by Sample\n")
        f.write("------------------------------------------\n")
        f.write("{:<30} {:<15} {:<15} {:<15}\n".format("Sample", "Total", "TC>TT", "GA>AA"))
        f.write("{:<30} {:<15} {:<15} {:<15}\n".format("-"*30, "-"*15, "-"*15, "-"*15))
        
        for sample, counts in sorted(sample_data.items()):
            f.write("{:<30} {:<15} {:<15} {:<15}\n".format(
                sample,
                counts['total'],
                counts['tc_tt'],
                counts['ga_aa']
            ))
        
        f.write("\n\n")
        
        # Calculate overall statistics
        total_mutations = sum(data['total'] for data in sample_data.values())
        tc_tt_mutations = sum(data['tc_tt'] for data in sample_data.values())
        ga_aa_mutations = sum(data['ga_aa'] for data in sample_data.values())
        
        f.write(f"Total APOBEC3-related mutations across all samples: {total_mutations}\n")
        f.write(f"Total TC>TT mutations: {tc_tt_mutations}\n")
        f.write(f"Total GA>AA mutations: {ga_aa_mutations}\n\n")
        
        # Per sample breakdown (detailed)
        f.write("Per Sample Breakdown (Detailed):\n")
        f.write("-------------------------------\n")
        for sample, counts in sorted(sample_data.items()):
            f.write(f"\nSample: {sample}\n")
            f.write(f"Total mutations: {counts['total']}\n")
            f.write(f"TC>TT mutations: {counts['tc_tt']}\n")
            f.write(f"GA>AA mutations: {counts['ga_aa']}\n")
    
    print(f"\nSummary statistics saved to apobec_summary.txt")

if __name__ == '__main__':
    main() 