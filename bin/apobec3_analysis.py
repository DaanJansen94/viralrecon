#!/usr/bin/env python

import pysam
import pandas as pd
import os
import glob
import sys

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
                sample_mutations.append(mutation_info_copy)
                all_apobec_mutations.append(mutation_info_copy)
                
        except Exception as e:
            print(f"Error processing variant at {chrom}:{pos}: {e}")

    # Create DataFrame for this sample and save to CSV
    sample_df = pd.DataFrame(sample_mutations)
    if not sample_df.empty:
        sample_output_file = f"{sample_name}_apobec_mutations.csv"
        
        # Calculate sample summary statistics
        tc_tt_count = len(sample_df[sample_df['type'] == 'TC>TT'])
        ga_aa_count = len(sample_df[sample_df['type'] == 'GA>AA'])
        
        # Save sample mutations to CSV
        sample_df.to_csv(sample_output_file, index=False)
        
        # Append summary statistics to the sample CSV file
        with open(sample_output_file, 'a') as f:
            f.write('\n\nSummary Statistics\n')
            f.write(f'Total APOBEC3-related mutations,{len(sample_mutations)}\n')
            f.write(f'TC>TT mutations,{tc_tt_count}\n')
            f.write(f'GA>AA mutations,{ga_aa_count}\n')
        
        print(f"Sample mutations saved to {sample_output_file}")
    
    # Print summary for this sample
    print(f"\nSummary for {vcf_file}:")
    print(f"Total TC>TT mutations found: {tc_tt_count}")
    print(f"Total GA>AA mutations found: {ga_aa_count}")
    print(f"Total APOBEC3-related mutations: {len(sample_mutations)}")

    # Load existing summary data if it exists
    summary_file = "apobec_summary.txt"
    existing_data = {}
    
    if os.path.exists(summary_file):
        try:
            # Try to extract existing sample data
            with open(summary_file, 'r') as f:
                lines = f.readlines()
                in_sample_section = False
                current_sample = None
                
                for line in lines:
                    if line.startswith("Sample:"):
                        current_sample = line.split("Sample:")[1].strip()
                        existing_data[current_sample] = {'total': 0, 'tc_tt': 0, 'ga_aa': 0}
                        in_sample_section = True
                    elif in_sample_section and line.startswith("Total mutations:"):
                        existing_data[current_sample]['total'] = int(line.split(":")[1].strip())
                    elif in_sample_section and line.startswith("TC>TT mutations:"):
                        existing_data[current_sample]['tc_tt'] = int(line.split(":")[1].strip())
                    elif in_sample_section and line.startswith("GA>AA mutations:"):
                        existing_data[current_sample]['ga_aa'] = int(line.split(":")[1].strip())
        except Exception as e:
            print(f"Error parsing existing summary file: {e}")
            # If there's an error parsing, we'll create a new file
            existing_data = {}

    # Process current sample mutations
    if sample_df.empty:
        return  # Skip if no mutations found
        
    # Current sample data
    sample_data = {
        'total': len(sample_mutations),
        'tc_tt': tc_tt_count,
        'ga_aa': ga_aa_count
    }
    
    # Add current sample to existing data
    existing_data[sample_name] = sample_data
    
    # Collect all sample mutations
    all_sample_files = glob.glob("*_apobec_mutations.csv")
    all_mutations = []
    
    for file in all_sample_files:
        try:
            # Extract sample name from filename
            file_sample = file.replace("_apobec_mutations.csv", "")
            
            # Only read the data part (not the summary at the end)
            df = pd.read_csv(file)
            if not df.empty:
                # If sample not in existing_data, add it
                if file_sample not in existing_data:
                    tc_tt = len(df[df['type'] == 'TC>TT'])
                    ga_aa = len(df[df['type'] == 'GA>AA'])
                    existing_data[file_sample] = {
                        'total': len(df),
                        'tc_tt': tc_tt,
                        'ga_aa': ga_aa
                    }
                all_mutations.append(df)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    # Combine all mutations if we have any
    if all_mutations:
        try:
            all_mutations_df = pd.concat(all_mutations, ignore_index=True)
        except Exception as e:
            print(f"Error combining mutation data: {e}")
            all_mutations_df = pd.DataFrame(all_apobec_mutations)
    else:
        all_mutations_df = pd.DataFrame(all_apobec_mutations)
    
    # Now write the updated summary file
    with open(summary_file, 'w') as f:
        f.write("APOBEC3 Mutation Analysis Summary\n")
        f.write("================================\n\n")
        
        # Create the summary table at the top
        f.write("Summary Table - APOBEC3 Mutations by Sample\n")
        f.write("------------------------------------------\n")
        f.write("{:<30} {:<15} {:<15} {:<15}\n".format("Sample", "Total", "TC>TT", "GA>AA"))
        f.write("{:<30} {:<15} {:<15} {:<15}\n".format("-"*30, "-"*15, "-"*15, "-"*15))
        
        for sample, counts in sorted(existing_data.items()):
            f.write("{:<30} {:<15} {:<15} {:<15}\n".format(
                sample,
                counts['total'],
                counts['tc_tt'],
                counts['ga_aa']
            ))
        
        f.write("\n\n")
        
        # Calculate overall statistics
        total_mutations = sum(data['total'] for data in existing_data.values())
        tc_tt_mutations = sum(data['tc_tt'] for data in existing_data.values())
        ga_aa_mutations = sum(data['ga_aa'] for data in existing_data.values())
        
        f.write(f"Total APOBEC3-related mutations across all samples: {total_mutations}\n")
        f.write(f"Total TC>TT mutations: {tc_tt_mutations}\n")
        f.write(f"Total GA>AA mutations: {ga_aa_mutations}\n\n")
        
        # Per sample breakdown (detailed)
        f.write("Per Sample Breakdown (Detailed):\n")
        f.write("-------------------------------\n")
        for sample, counts in sorted(existing_data.items()):
            f.write(f"\nSample: {sample}\n")
            f.write(f"Total mutations: {counts['total']}\n")
            f.write(f"TC>TT mutations: {counts['tc_tt']}\n")
            f.write(f"GA>AA mutations: {counts['ga_aa']}\n")
        
    print(f"\nUpdated summary statistics saved to {summary_file}")

if __name__ == '__main__':
    main() 