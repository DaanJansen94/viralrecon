#!/usr/bin/env python

import gzip
import sys
import os
import csv

def count_reads_in_fastq(fastq_file):
    """Count the number of reads in a FASTQ file (gzipped or regular)"""
    line_count = 0
    
    # Check if file exists
    if not os.path.exists(fastq_file):
        print(f"Error: File {fastq_file} not found")
        return 0
    
    # Open file (gzipped or regular)
    opener = gzip.open if fastq_file.endswith('.gz') else open
    try:
        with opener(fastq_file, 'rt') as f:
            for line in f:
                line_count += 1
    except Exception as e:
        print(f"Error reading {fastq_file}: {e}")
        return 0
    
    # Each read in FASTQ has 4 lines
    return line_count // 4

def main():
    """Main function to parse samplesheet and count reads"""
    if len(sys.argv) != 2:
        print("Usage: count_fastq_reads.py <samplesheet.csv>")
        sys.exit(1)
    
    samplesheet = sys.argv[1]
    output_file = "read_counts.csv"
    
    # Check if samplesheet exists
    if not os.path.exists(samplesheet):
        print(f"Error: Samplesheet {samplesheet} not found")
        sys.exit(1)
    
    # Read the samplesheet and count reads
    sample_read_counts = {}
    
    try:
        with open(samplesheet, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                sample = row.get('sample', '')
                fastq_1 = row.get('fastq_1', '')
                fastq_2 = row.get('fastq_2', '')
                
                if not (sample and fastq_1):
                    print(f"Warning: Invalid row in samplesheet: {row}")
                    continue
                
                # Count reads in R1 file (for paired-end, R1 and R2 should have the same number)
                read_count = count_reads_in_fastq(fastq_1)
                
                # For paired-end data, the total read count is the number of pairs
                sample_read_counts[sample] = read_count
                
                print(f"Sample: {sample}, Read count: {read_count}")
    
    except Exception as e:
        print(f"Error processing samplesheet: {e}")
        sys.exit(1)
    
    # Write results to output file
    try:
        with open(output_file, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['sample', 'raw_read_count'])
            for sample, count in sample_read_counts.items():
                writer.writerow([sample, count])
        
        print(f"Read counts written to {output_file}")
    
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 