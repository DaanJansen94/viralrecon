#!/usr/bin/env python

import os
import glob
import sys
import csv

def main():
    if len(sys.argv) != 2:
        print("Usage: combine_apobec_summaries.py <output_dir>")
        sys.exit(1)
    
    output_dir = sys.argv[1]
    print(f"Processing APOBEC3 summaries in: {output_dir}")
    
    # Find all apobec CSV files directly
    apobec_csv_files = glob.glob(os.path.join(output_dir, "*_apobec_mutations.csv"))
    
    print(f"Found {len(apobec_csv_files)} APOBEC3 result files")
    
    # Extract data from each CSV file
    sample_data = {}
    
    for file_path in apobec_csv_files:
        try:
            sample_name = os.path.basename(file_path).replace("_apobec_mutations.csv", "")
            
            # Initialize with default values
            data = {'all': 0, 'total': 0, 'tc_tt': 0, 'ga_aa': 0}
            
            # Read the summary statistics section from the CSV file
            with open(file_path, 'r') as f:
                lines = f.readlines()
                read_count = None
                
                for i, line in enumerate(lines):
                    if 'Total mutations,' in line:
                        try:
                            data['all'] = int(line.strip().split(',')[1])
                        except:
                            print(f"Warning: Could not parse total mutations from {file_path}")
                    elif 'Total APOBEC3-related mutations,' in line:
                        try:
                            data['total'] = int(line.strip().split(',')[1])
                        except:
                            print(f"Warning: Could not parse APOBEC3 mutations from {file_path}")
                    elif 'TC>TT mutations,' in line:
                        try:
                            data['tc_tt'] = int(line.strip().split(',')[1])
                        except:
                            print(f"Warning: Could not parse TC>TT mutations from {file_path}")
                    elif 'GA>AA mutations,' in line:
                        try:
                            data['ga_aa'] = int(line.strip().split(',')[1])
                        except:
                            print(f"Warning: Could not parse GA>AA mutations from {file_path}")
                    elif 'Raw read count,' in line:
                        try:
                            read_count = int(line.strip().split(',')[1])
                            data['read_count'] = read_count
                        except:
                            print(f"Warning: Could not parse raw read count from {file_path}")
            
            # Add the sample data
            sample_data[sample_name] = data
            print(f"Processed {sample_name}: {data['all']} mutations, {data['total']} APOBEC3 mutations" + 
                  (f", {read_count} reads" if read_count else ""))
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
    
    # Also load read counts from separate file if available, to catch any that were missed
    read_counts_file = os.path.join(output_dir, "read_counts.csv")
    if os.path.exists(read_counts_file):
        try:
            with open(read_counts_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if 'sample' in row and 'raw_read_count' in row:
                        sample_name = row['sample']
                        read_count = int(row['raw_read_count'])
                        
                        # If we already have this sample, add the read count
                        if sample_name in sample_data and 'read_count' not in sample_data[sample_name]:
                            sample_data[sample_name]['read_count'] = read_count
                            print(f"Added read count to existing sample {sample_name}: {read_count}")
                        # If we don't have this sample yet, add it with just the read count
                        elif sample_name not in sample_data:
                            sample_data[sample_name] = {
                                'all': 0,
                                'total': 0,
                                'tc_tt': 0,
                                'ga_aa': 0,
                                'read_count': read_count
                            }
                            print(f"Added sample with only read count: {sample_name}: {read_count}")
        except Exception as e:
            print(f"Error reading read counts file: {e}")
    else:
        print(f"Warning: Read counts file not found at {read_counts_file}")
    
    # Check if we have any samples at all
    if not sample_data:
        print("Warning: No sample data found. Creating minimal summary file.")
        # Create a minimal summary file for at least one sample
        read_counts = {}
        try:
            with open(read_counts_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if 'sample' in row and 'raw_read_count' in row:
                        sample_name = row['sample']
                        read_count = int(row['raw_read_count'])
                        read_counts[sample_name] = read_count
                        
            # Add at least the first sample from read counts
            if read_counts:
                sample_name = next(iter(read_counts))
                sample_data[sample_name] = {
                    'all': 0,
                    'total': 0,
                    'tc_tt': 0,
                    'ga_aa': 0,
                    'read_count': read_counts[sample_name]
                }
        except:
            print("Error creating minimal sample data")
    
    # Determine if we have read counts for any samples
    has_read_counts = any('read_count' in data for data in sample_data.values())
    
    # Create the summary file
    summary_file = os.path.join(output_dir, "apobec_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("APOBEC3 Mutation Analysis Summary\n")
        f.write("================================\n\n")
        
        # Add note about SNPs vs indels at the top
        f.write("NOTE: This analysis focuses on single nucleotide polymorphisms (SNPs) only.\n")
        f.write("Results do not include insertions, deletions, or other complex variants.\n")
        f.write("'All Mutations' refers to the total number of SNPs detected.\n")
        f.write("'APM' stands for APOBEC3 mutations Per Million reads.\n\n")
        
        # Create the summary table at the top
        f.write("Summary Table - APOBEC3 Mutations by Sample\n")
        f.write("------------------------------------------\n")
        if has_read_counts:
            f.write("{:<30} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n".format(
                "Sample", "Raw Reads", "All Mutations", "APOBEC3", "APOBEC3 (%)", "APM", "TC>TT", "GA>AA"))
            f.write("{:<30} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n".format(
                "-"*30, "-"*15, "-"*15, "-"*15, "-"*15, "-"*15, "-"*15, "-"*15))
        else:
            f.write("{:<30} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n".format(
                "Sample", "All Mutations", "APOBEC3", "APOBEC3 (%)", "APM", "TC>TT", "GA>AA"))
            f.write("{:<30} {:<15} {:<15} {:<15} {:<15} {:<15} {:<15}\n".format(
                "-"*30, "-"*15, "-"*15, "-"*15, "-"*15, "-"*15, "-"*15))
        
        for sample, counts in sorted(sample_data.items()):
            # Calculate percentage of APOBEC3 mutations
            all_muts = counts.get('all', 0)
            total_apobec = counts.get('total', 0)
            apobec_pct = 0
            if all_muts > 0:
                apobec_pct = (total_apobec / all_muts) * 100
            
            # Calculate APOBEC3 mutations per million reads (APM)
            apm = 0
            if 'read_count' in counts and counts['read_count'] > 0:
                apm = (total_apobec / counts['read_count']) * 1000000
            
            if has_read_counts:
                f.write("{:<30} {:<15} {:<15} {:<15} {:<15.2f} {:<15.2f} {:<15} {:<15}\n".format(
                    sample,
                    counts.get('read_count', 'N/A'),
                    all_muts,
                    total_apobec,
                    apobec_pct,
                    apm,
                    counts.get('tc_tt', 0),
                    counts.get('ga_aa', 0)
                ))
            else:
                f.write("{:<30} {:<15} {:<15} {:<15.2f} {:<15.2f} {:<15} {:<15}\n".format(
                    sample,
                    all_muts,
                    total_apobec,
                    apobec_pct,
                    apm,
                    counts.get('tc_tt', 0),
                    counts.get('ga_aa', 0)
                ))
        
        f.write("\n")
        
        # Calculate overall statistics
        total_samples = len(sample_data)
        total_all_mutations = sum(data.get('all', 0) for data in sample_data.values())
        total_mutations = sum(data.get('total', 0) for data in sample_data.values())
        total_tc_tt = sum(data.get('tc_tt', 0) for data in sample_data.values())
        total_ga_aa = sum(data.get('ga_aa', 0) for data in sample_data.values())
        
        # Calculate total percentage
        overall_pct = 0
        if total_all_mutations > 0:
            overall_pct = (total_mutations / total_all_mutations) * 100
        
        # Calculate overall APM
        overall_apm = 0
        if has_read_counts:
            total_reads = sum(data.get('read_count', 0) for data in sample_data.values() if 'read_count' in data)
            if total_reads > 0:
                overall_apm = (total_mutations / total_reads) * 1000000
        
        f.write(f"\nAnalysis Summary:\n")
        f.write(f"Total samples analyzed: {total_samples}\n")
        f.write(f"Total mutations (all samples): {total_all_mutations}\n")
        if has_read_counts:
            total_reads = sum(data.get('read_count', 0) for data in sample_data.values() if 'read_count' in data)
            f.write(f"Total raw reads (all samples): {total_reads}\n")
        f.write(f"Total APOBEC3-related mutations: {total_mutations} ({overall_pct:.1f}% of all mutations)\n")
        if has_read_counts and total_reads > 0:
            f.write(f"Overall APOBEC3 mutations per million reads (APM): {overall_apm:.2f}\n")
        f.write(f"Total TC>TT mutations: {total_tc_tt}\n")
        f.write(f"Total GA>AA mutations: {total_ga_aa}\n\n")
        
        # Per-sample detailed breakdown
        f.write("Individual Sample Details:\n")
        f.write("------------------------\n")
        
        for sample, counts in sorted(sample_data.items()):
            f.write(f"\nSample: {sample}\n")
            if 'read_count' in counts:
                f.write(f"Raw reads: {counts['read_count']}\n")
            
            all_muts = counts.get('all', 0)
            total_apobec = counts.get('total', 0)
            apobec_pct = 0
            if all_muts > 0:
                apobec_pct = (total_apobec / all_muts) * 100
            
            # Calculate APOBEC3 mutations per million reads
            apm = 0
            if 'read_count' in counts and counts['read_count'] > 0:
                apm = (total_apobec / counts['read_count']) * 1000000
            
            f.write(f"All mutations: {all_muts}\n")
            f.write(f"APOBEC3 mutations: {total_apobec} ({apobec_pct:.1f}% of all mutations)\n")
            if 'read_count' in counts:
                f.write(f"APOBEC3 mutations per million reads (APM): {apm:.2f}\n")
            f.write(f"TC>TT mutations: {counts.get('tc_tt', 0)}\n")
            f.write(f"GA>AA mutations: {counts.get('ga_aa', 0)}\n")
    
    print(f"Summary written to: {summary_file}")

if __name__ == "__main__":
    main() 