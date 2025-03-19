#!/usr/bin/env python

import os
import sys
import glob

def update_summary_file(file_path):
    # Read the current content
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Check if the note is already present to avoid duplication
    if 'Note: This analysis specifically focuses on single nucleotide polymorphisms' not in content:
        # Split at the first heading after the title
        parts = content.split('Summary Table', 1)
        if len(parts) == 2:
            # Add the note before the summary table
            new_content = parts[0]
            new_content += "Note: This analysis specifically focuses on single nucleotide polymorphisms (SNPs)\n"
            new_content += "and does not include insertions, deletions, or other complex variants.\n"
            new_content += "'All Mutations' refers to the total number of SNPs detected.\n\n"
            new_content += "Summary Table" + parts[1]
            
            # Write the updated content back
            with open(file_path, 'w') as f:
                f.write(new_content)
            print(f"Updated {file_path} with SNP note")
            return True
    
    print(f"Note already exists or file format not as expected in {file_path}")
    return False

def main():
    # Get the apobec summary file path
    if len(sys.argv) > 1:
        summary_file = sys.argv[1]
    else:
        # Try to find the file in the current directory
        files = glob.glob("**/apobec_summary.txt", recursive=True)
        if files:
            summary_file = files[0]
        else:
            print("No apobec_summary.txt file found. Please provide path as argument.")
            return
    
    # Update the file
    if os.path.exists(summary_file):
        update_summary_file(summary_file)
    else:
        print(f"File not found: {summary_file}")

if __name__ == "__main__":
    main() 