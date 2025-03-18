#!/usr/bin/env python

import os
import glob
import sys

def main():
    if len(sys.argv) != 2:
        print("Usage: combine_apobec_summaries.py <output_dir>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    print(f"Processing APOBEC3 summaries in: {output_dir}")
    
    # Find all sample summary data files
    sample_data_files = glob.glob(os.path.join(output_dir, "*_summary_data.txt"))
    apobec_csv_files = glob.glob(os.path.join(output_dir, "*_apobec_mutations.csv"))
    
    print(f"Found {len(sample_data_files)} sample summary files")
    print(f"Found {len(apobec_csv_files)} APOBEC3 result files")
    
    # Extract data from each summary file
    sample_data = {}
    
    for file_path in sample_data_files:
        try:
            sample_name = None
            data = {'total': 0, 'tc_tt': 0, 'ga_aa': 0}
            
            with open(file_path, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith("Sample:"):
                        sample_name = line.split("Sample:")[1].strip()
                    elif line.startswith("Total:"):
                        data['total'] = int(line.split("Total:")[1].strip())
                    elif line.startswith("TC>TT:"):
                        data['tc_tt'] = int(line.split("TC>TT:")[1].strip())
                    elif line.startswith("GA>AA:"):
                        data['ga_aa'] = int(line.split("GA>AA:")[1].strip())
            
            if sample_name:
                sample_data[sample_name] = data
                print(f"Loaded summary data for {sample_name}")
            else:
                print(f"Warning: Could not extract sample name from {file_path}")
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    # Create the combined summary file
    summary_file = os.path.join(output_dir, "apobec_combined_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("APOBEC3 Mutation Analysis Combined Summary\n")
        f.write("========================================\n\n")
        
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
    
    print(f"Combined summary written to: {summary_file}")
    
    # Copy the combined summary to also replace the regular summary file
    regular_summary = os.path.join(output_dir, "apobec_summary.txt")
    try:
        with open(summary_file, 'r') as src:
            with open(regular_summary, 'w') as dest:
                dest.write(src.read())
        print(f"Also copied to standard summary file: {regular_summary}")
    except Exception as e:
        print(f"Error copying to standard summary file: {e}")

if __name__ == "__main__":
    main() 