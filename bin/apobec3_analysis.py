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
        
        # Initialize empty variables
        tc_tt_count = 0
        ga_aa_count = 0
        total_mutations = 0
        total_all_mutations = 0
        
        # Create empty summary files
        sample_output_file = f"{sample_name}_apobec_mutations.csv"
        with open(sample_output_file, 'w') as f:
            f.write('chrom,pos,ref,alt,depth,quality,sample,context,type\n\n')
            f.write('Summary Statistics\n')
            f.write(f'Total mutations,{total_all_mutations}\n')
            f.write(f'Total APOBEC3-related mutations,{total_mutations}\n')
            f.write(f'TC>TT mutations,{tc_tt_count}\n')
            f.write(f'GA>AA mutations,{ga_aa_count}\n')
        
        # Write simple summary file
        with open("apobec_summary.txt", 'w') as f:
            f.write("APOBEC3 Mutation Analysis Summary\n")
            f.write("================================\n\n")
            f.write("Note: This analysis specifically focuses on single nucleotide polymorphisms (SNPs)\n")
            f.write("and does not include insertions, deletions, or other complex variants.\n")
            f.write("'All Mutations' refers to the total number of SNPs detected.\n\n")
            f.write(f"Error processing {vcf_file}: {e}\n\n")
            f.write("No mutations were analyzed due to file access error.\n")
        
        print(f"Created empty summary files due to error: {e}")
        sys.exit(1)

    # Initialize storage for this sample's mutations
    sample_mutations = []
    tc_tt_mutations = []
    ga_aa_mutations = []

    # Count total mutations in the VCF file
    total_all_mutations = 0
    for variant in vcf:
        if variant.alts and len(variant.ref) == 1 and len(variant.alts[0]) == 1:  # Count only SNPs
            total_all_mutations += 1
    
    # Reset the VCF reader to start from the beginning
    vcf = pysam.VariantFile(vcf_file)

    # Initialize count variables
    tc_tt_count = 0
    ga_aa_count = 0
    total_mutations = 0

    # Iterate through variants in the VCF file to find APOBEC3 mutations
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
    
    # Calculate sample summary statistics
    tc_tt_count = len(tc_tt_mutations)
    ga_aa_count = len(ga_aa_mutations)
    total_mutations = tc_tt_count + ga_aa_count
    
    if not sample_df.empty:
        sample_output_file = f"{sample_name}_apobec_mutations.csv"
        
        # Save sample mutations to CSV
        sample_df.to_csv(sample_output_file, index=False)
        
        # Append summary statistics to the sample CSV file
        with open(sample_output_file, 'a') as f:
            f.write('\n\nSummary Statistics\n')
            f.write(f'Total mutations,{total_all_mutations}\n')
            f.write(f'Total APOBEC3-related mutations,{total_mutations}\n')
            f.write(f'TC>TT mutations,{tc_tt_count}\n')
            f.write(f'GA>AA mutations,{ga_aa_count}\n')
        
        print(f"Sample mutations saved to {sample_output_file}")
    else:
        # Create empty file with just summary statistics
        sample_output_file = f"{sample_name}_apobec_mutations.csv"
        with open(sample_output_file, 'w') as f:
            f.write('chrom,pos,ref,alt,depth,quality,sample,context,type\n\n')
            f.write('Summary Statistics\n')
            f.write(f'Total mutations,{total_all_mutations}\n')
            f.write(f'Total APOBEC3-related mutations,{total_mutations}\n')
            f.write(f'TC>TT mutations,{tc_tt_count}\n')
            f.write(f'GA>AA mutations,{ga_aa_count}\n')
        print(f"No APOBEC3 mutations found. Empty CSV with summary created: {sample_output_file}")
    
    # Print summary for this sample
    print(f"\nSummary for {vcf_file}:")
    print(f"Total SNPs: {total_all_mutations}")
    print(f"Total TC>TT mutations found: {tc_tt_count}")
    print(f"Total GA>AA mutations found: {ga_aa_count}")
    print(f"Total APOBEC3-related mutations: {total_mutations}")

    # Create a summary file
    sample_data = {}
    
    # Process any existing result files and add our own
    current_files = glob.glob("*_apobec_mutations.csv")
    print(f"Found {len(current_files)} sample file(s) in current directory")
    
    for file in current_files:
        try:
            file_sample = os.path.basename(file).replace("_apobec_mutations.csv", "")
            
            # Try to read summary data from file
            with open(file, 'r') as f:
                lines = f.readlines()
                all_muts = 0
                apobec_muts = 0
                tc_tt_muts = 0
                ga_aa_muts = 0
                
                for i, line in enumerate(lines):
                    if 'Total mutations,' in line:
                        all_muts = int(line.strip().split(',')[1])
                    elif 'Total APOBEC3-related mutations,' in line:
                        apobec_muts = int(line.strip().split(',')[1])
                    elif 'TC>TT mutations,' in line:
                        tc_tt_muts = int(line.strip().split(',')[1])
                    elif 'GA>AA mutations,' in line:
                        ga_aa_muts = int(line.strip().split(',')[1])
                
                sample_data[file_sample] = {
                    'all': all_muts,
                    'total': apobec_muts, 
                    'tc_tt': tc_tt_muts, 
                    'ga_aa': ga_aa_muts
                }
                
                print(f"Added sample {file_sample} with {all_muts} total mutations, {apobec_muts} APOBEC3 mutations")
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    # Write our current sample data if not already in the data
    if sample_name not in sample_data:
        sample_data[sample_name] = {
            'all': total_all_mutations,
            'total': total_mutations,
            'tc_tt': tc_tt_count,
            'ga_aa': ga_aa_count
        }
    
    # Write summary file
    with open("apobec_summary.txt", 'w') as f:
        f.write("APOBEC3 Mutation Analysis Summary\n")
        f.write("================================\n\n")
        
        # Add clarification about analysis
        f.write("Note: This analysis specifically focuses on single nucleotide polymorphisms (SNPs)\n")
        f.write("and does not include insertions, deletions, or other complex variants.\n")
        f.write("'All Mutations' refers to the total number of SNPs detected.\n\n")
        
        # Create the summary table
        f.write("Summary Table - APOBEC3 Mutations by Sample\n")
        f.write("------------------------------------------\n")
        f.write("{:<30} {:<15} {:<20} {:<15} {:<15}\n".format("Sample", "All Mutations", "APOBEC3 (%)", "TC>TT", "GA>AA"))
        f.write("{:<30} {:<15} {:<20} {:<15} {:<15}\n".format("-"*30, "-"*15, "-"*20, "-"*15, "-"*15))
        
        for sample, counts in sorted(sample_data.items()):
            # Calculate percentage of APOBEC3 mutations
            all_muts = counts.get('all', 0)
            apobec_pct = 0
            if all_muts > 0:
                apobec_pct = (counts['total'] / all_muts) * 100
            
            f.write("{:<30} {:<15} {:<20} {:<15} {:<15}\n".format(
                sample,
                counts.get('all', 0),
                f"{counts['total']} ({apobec_pct:.1f}%)",
                counts['tc_tt'],
                counts['ga_aa']
            ))
        
        f.write("\n\n")
        
        # Calculate overall statistics
        total_all_muts = sum(data.get('all', 0) for data in sample_data.values())
        total_mutations = sum(data['total'] for data in sample_data.values())
        tc_tt_mutations = sum(data['tc_tt'] for data in sample_data.values())
        ga_aa_mutations = sum(data['ga_aa'] for data in sample_data.values())
        
        # Calculate overall percentage
        overall_pct = 0
        if total_all_muts > 0:
            overall_pct = (total_mutations / total_all_muts) * 100
        
        f.write(f"Total mutations across all samples: {total_all_muts}\n")
        f.write(f"Total APOBEC3-related mutations: {total_mutations} ({overall_pct:.1f}% of all mutations)\n")
        f.write(f"Total TC>TT mutations: {tc_tt_mutations}\n")
        f.write(f"Total GA>AA mutations: {ga_aa_mutations}\n\n")
        
        # Per sample breakdown (detailed)
        f.write("Per Sample Breakdown (Detailed):\n")
        f.write("-------------------------------\n")
        for sample, counts in sorted(sample_data.items()):
            # Calculate percentage for individual sample
            all_muts = counts.get('all', 0)
            apobec_pct = 0
            if all_muts > 0:
                apobec_pct = (counts['total'] / all_muts) * 100
            
            f.write(f"\nSample: {sample}\n")
            f.write(f"All mutations: {counts.get('all', 0)}\n")
            f.write(f"APOBEC3 mutations: {counts['total']} ({apobec_pct:.1f}% of all mutations)\n")
            f.write(f"TC>TT mutations: {counts['tc_tt']}\n")
            f.write(f"GA>AA mutations: {counts['ga_aa']}\n")
    
    print(f"\nSummary statistics saved to apobec_summary.txt")

if __name__ == '__main__':
    main() 