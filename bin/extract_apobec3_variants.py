#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import pandas as pd
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Extract APOBEC3 variants (G>A and C>T mutations) from iVar TSV files')
    parser.add_argument('--input_dir', required=True, help='Directory containing iVar variant TSV files')
    parser.add_argument('--output_dir', required=True, help='Directory to output filtered APOBEC3 variant files')
    parser.add_argument('--min_freq', type=float, default=0.0, help='Minimum allele frequency threshold (default: 0.0)')
    parser.add_argument('--min_depth', type=int, default=0, help='Minimum depth threshold (default: 0)')
    parser.add_argument('--major_cutoff', type=float, default=0.5, help='Cutoff for major variants (default: 0.5)')
    return parser.parse_args()

def is_apobec3_mutation(row):
    """
    Determine if a variant is an APOBEC3-type mutation (G>A or C>T)
    """
    return ((row['REF'] == 'G' and row['ALT'] == 'A') or 
            (row['REF'] == 'C' and row['ALT'] == 'T'))

def generate_frequency_distribution(df):
    """
    Generate a table showing the distribution of variants across frequency ranges
    """
    # Define frequency bins
    bins = [0, 0.005, 0.01, 0.05, 0.1, 0.5]
    bin_labels = ['0-0.5%', '0.5-1%', '1-5%', '5-10%', '10-50%']
    
    # Count variants in each bin
    freq_counts = pd.cut(df['ALT_FREQ'], bins=bins, labels=bin_labels).value_counts().sort_index()
    
    # Create a DataFrame for the summary
    summary_df = pd.DataFrame({
        'Frequency_Range': freq_counts.index,
        'Variant_Count': freq_counts.values
    })
    
    return summary_df

def process_file(input_file, output_base, min_freq=0.0, min_depth=0, major_cutoff=0.5):
    """
    Process an iVar TSV file and extract APOBEC3 variants,
    separating them into major and minor variant files
    """
    try:
        # Read the TSV file
        df = pd.read_csv(input_file, sep='\t')
        
        # Extract sample name from file path
        sample_name = os.path.basename(output_base)
        
        # Filter for APOBEC3 mutations (G>A or C>T)
        apobec3_df = df[df.apply(is_apobec3_mutation, axis=1)]
        
        # Apply frequency and depth filters
        if min_freq > 0:
            apobec3_df = apobec3_df[apobec3_df['ALT_FREQ'] >= min_freq]
        
        if min_depth > 0:
            apobec3_df = apobec3_df[apobec3_df['TOTAL_DP'] >= min_depth]
        
        # Split into major and minor variants
        major_variants = apobec3_df[apobec3_df['ALT_FREQ'] >= major_cutoff]
        minor_variants = apobec3_df[apobec3_df['ALT_FREQ'] < major_cutoff]
        
        # Generate frequency distribution for minor variants
        freq_dist = None
        if not minor_variants.empty:
            freq_dist = generate_frequency_distribution(minor_variants)
            
            # Create a string representation of the frequency distribution table
            freq_table = "# Frequency Distribution of Minor APOBEC3 Variants\n"
            freq_table += "# " + "-" * 50 + "\n"
            freq_table += "# Frequency_Range\tVariant_Count\n"
            
            for _, row in freq_dist.iterrows():
                freq_table += f"# {row['Frequency_Range']}\t{row['Variant_Count']}\n"
            
            freq_table += "# " + "-" * 50 + "\n\n"
        else:
            freq_table = "# No minor APOBEC3 variants found\n\n"
        
        # Create output filenames
        all_output = f"{output_base}.apobec3.tsv"
        major_output = f"{output_base}.apobec3.major.tsv"
        minor_output = f"{output_base}.apobec3.minor.tsv"
        
        # Write filtered data to output files
        apobec3_df.to_csv(all_output, sep='\t', index=False)
        major_variants.to_csv(major_output, sep='\t', index=False)
        
        # For minor variants, prepend the frequency distribution table
        with open(minor_output, 'w') as f:
            f.write(freq_table)
            if not minor_variants.empty:
                minor_variants.to_csv(f, sep='\t', index=False)
            else:
                f.write("# No data\n")
        
        print(f"Processed {input_file}:")
        print(f"  Total variants: {len(df)}")
        print(f"  APOBEC3 variants: {len(apobec3_df)}")
        print(f"  Major APOBEC3 variants (≥{major_cutoff}): {len(major_variants)}")
        print(f"  Minor APOBEC3 variants (<{major_cutoff}): {len(minor_variants)}")
        
        # Print frequency distribution
        if not minor_variants.empty:
            print("\nFrequency Distribution of Minor APOBEC3 Variants:")
            for _, row in freq_dist.iterrows():
                print(f"  {row['Frequency_Range']}: {row['Variant_Count']} variants")
        
        # Return the frequency distribution with sample name for the summary
        if freq_dist is not None:
            freq_dist['Sample'] = sample_name
            return len(apobec3_df), freq_dist
        else:
            # Create an empty dataframe with the same structure
            empty_dist = pd.DataFrame({
                'Frequency_Range': bin_labels,
                'Variant_Count': [0] * len(bin_labels),
                'Sample': [sample_name] * len(bin_labels)
            })
            return len(apobec3_df), empty_dist
    except Exception as e:
        print(f"Error processing {input_file}: {str(e)}")
        # Return empty dataframe with the right structure for summary
        empty_dist = pd.DataFrame({
            'Frequency_Range': bin_labels,
            'Variant_Count': [0] * len(bin_labels),
            'Sample': [os.path.basename(output_base)] * len(bin_labels)
        })
        return 0, empty_dist

def create_summary_file(all_freq_dists, output_dir):
    """
    Create a summary file with a table showing frequency distributions across all samples
    Samples will be displayed as rows and frequency ranges as columns
    """
    if not all_freq_dists:
        print("No frequency distributions available to create summary")
        return
    
    # Combine all frequency distributions
    combined_df = pd.concat(all_freq_dists)
    
    # Define the bin labels to ensure consistent order
    bin_labels = ['0-0.5%', '0.5-1%', '1-5%', '5-10%', '10-50%']
    
    # Create a DataFrame to store the results (samples as rows, frequency ranges as columns)
    samples = combined_df['Sample'].unique()
    result_df = pd.DataFrame(columns=['Sample'] + bin_labels)
    
    # Create rows for each sample
    for i, sample in enumerate(sorted(samples)):
        sample_data = combined_df[combined_df['Sample'] == sample]
        row_data = {'Sample': sample}
        
        # Get counts for each frequency range
        for freq_range in bin_labels:
            row = sample_data[sample_data['Frequency_Range'] == freq_range]
            if not row.empty:
                row_data[freq_range] = row.iloc[0]['Variant_Count']
            else:
                row_data[freq_range] = 0
                
        # Add the row to the result DataFrame
        result_df.loc[i] = row_data
    
    # Add a total column
    result_df['Total'] = result_df[bin_labels].sum(axis=1)
    
    # Create output file
    output_file = os.path.join(output_dir, "combined_summary_minor_variants.tsv")
    
    # Format the table to be more readable
    # Open the file and write a header
    with open(output_file, 'w') as f:
        # Write header
        f.write("# Summary of Minor APOBEC3 Variants Across Samples\n")
        f.write("# Each cell shows the number of variants in the given frequency range for each sample\n")
        f.write("#\n")
        
        # Get column widths for proper alignment
        col_widths = {col: max(len(str(col)), max(len(str(val)) for val in result_df[col])) + 2 
                     for col in result_df.columns}
        
        # Write header row
        for col in result_df.columns:
            f.write(f"{col:<{col_widths[col]}}")
        f.write("\n")
        
        # Write data rows
        for _, row in result_df.iterrows():
            for col in result_df.columns:
                f.write(f"{row[col]:<{col_widths[col]}}")
            f.write("\n")
    
    sample_count = len(samples)
    print(f"\nSummary of minor APOBEC3 variants by frequency range across {sample_count} samples:")
    print(f"Saved to: {output_file}")
    
    return output_file

def main():
    # Parse command line arguments
    args = parse_args()
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Get list of TSV files
    tsv_files = glob.glob(os.path.join(args.input_dir, "*.tsv"))
    
    if not tsv_files:
        print(f"No TSV files found in {args.input_dir}")
        sys.exit(1)
    
    print(f"Found {len(tsv_files)} TSV files")
    
    # Process each file
    total_apobec3 = 0
    all_freq_dists = []
    
    for input_file in tsv_files:
        filename = os.path.basename(input_file)
        output_base = os.path.join(args.output_dir, os.path.splitext(filename)[0])
        variants_count, freq_dist = process_file(
            input_file, 
            output_base, 
            args.min_freq, 
            args.min_depth,
            args.major_cutoff
        )
        total_apobec3 += variants_count
        all_freq_dists.append(freq_dist)
    
    # Create summary file with frequency distributions across all samples
    create_summary_file(all_freq_dists, args.output_dir)
    
    print(f"Total APOBEC3 variants found: {total_apobec3}")
    print(f"Results written to {args.output_dir}")

if __name__ == "__main__":
    main() 