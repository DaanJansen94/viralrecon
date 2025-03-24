process COMBINE_APOBEC3_VARIANTS {
    tag "combine_samples"
    label 'process_low'

    conda "conda-forge::python=3.9.5 conda-forge::pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9.5' :
        'python:3.9.5' }"

    publishDir "${params.outdir}/variants/ivar/apobec3", mode: params.publish_dir_mode, overwrite: true

    input:
    path(minor_variants)

    output:
    path "combined_summary_minor_variants.tsv", emit: summary_tsv
    path "combined_summary_minor_variants.txt", emit: summary_txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Set up temporary directory for pip packages
    mkdir -p /tmp/package_install
    export PIP_TARGET=/tmp/package_install
    export PYTHONPATH="/tmp/package_install"
    
    # Install required packages
    pip install --no-cache-dir --target=/tmp/package_install pandas==1.3.5 numpy==1.21.6
    
    # Create a script to combine the minor variant files
    cat > combine_summary.py << 'EOF'
#!/usr/bin/env python3
import os
import sys
import glob
import re
import pandas as pd

# Define frequency bins for consistency
bins = ['0-0.5%', '0.5-1%', '1-5%', '5-10%', '10-50%']

# Dictionary to store counts for each sample
sample_counts = {}

# Process each minor variant file
print("Looking for minor variant files...")
minor_files = glob.glob("*.apobec3.minor.tsv")
if not minor_files:
    print("No minor variant files found")
    sys.exit(1)
    
print(f"Found {len(minor_files)} minor variant files:")

for file in minor_files:
    # Extract sample name from filename
    filename = os.path.basename(file)
    sample_match = re.search(r'(.+?)\\.apobec3\\.minor\\.tsv', filename)
    if sample_match:
        sample_name = sample_match.group(1)
    else:
        sample_name = filename.replace(".apobec3.minor.tsv", "")
    
    print(f"  - Processing {sample_name} from {filename}")
    
    # Read the file, skip the header comments
    with open(file, 'r') as f:
        lines = f.readlines()
        data_start = 0
        for i, line in enumerate(lines):
            if not line.startswith('#'):
                data_start = i
                break
    
    # Always read the data directly for most accurate counting
    try:
        df = pd.read_csv(file, sep='\\t', skiprows=data_start)
        if 'ALT_FREQ' in df.columns:
            # Print some diagnostics
            print(f"    Read {len(df)} rows, ALT_FREQ range: {df['ALT_FREQ'].min():.6f} - {df['ALT_FREQ'].max():.6f}")
            
            # Count variants in each bin
            counts = {
                '0-0.5%': len(df[df['ALT_FREQ'] < 0.005]),
                '0.5-1%': len(df[(df['ALT_FREQ'] >= 0.005) & (df['ALT_FREQ'] < 0.01)]),
                '1-5%': len(df[(df['ALT_FREQ'] >= 0.01) & (df['ALT_FREQ'] < 0.05)]),
                '5-10%': len(df[(df['ALT_FREQ'] >= 0.05) & (df['ALT_FREQ'] < 0.1)]),
                '10-50%': len(df[(df['ALT_FREQ'] >= 0.1) & (df['ALT_FREQ'] <= 0.5)])
            }
            
            # Print detailed counts for debugging
            count_10_50 = len(df[(df['ALT_FREQ'] >= 0.1) & (df['ALT_FREQ'] <= 0.5)])
            freqs_10_50 = df.loc[(df['ALT_FREQ'] >= 0.1) & (df['ALT_FREQ'] <= 0.5), 'ALT_FREQ'].tolist()
            print(f"    Found {count_10_50} variants in 10-50% range: {freqs_10_50}")
            
            sample_counts[sample_name] = counts
            print(f"    Found {sum(counts.values())} minor variants total")
        else:
            print(f"    Warning: ALT_FREQ column not found in {file}")
    except Exception as e:
        print(f"    Error reading variant data: {str(e)}")
        # Default to empty counts if we can't read the data
        sample_counts[sample_name] = {bin_name: 0 for bin_name in bins}

# Create the output DataFrame with samples as rows and frequency ranges as columns
result_df = pd.DataFrame(columns=['Sample'] + bins)

# Sort sample names for consistent output
sorted_samples = sorted(sample_counts.keys())

# Add each sample as a row
row_idx = 0
for sample in sorted_samples:
    counts = sample_counts[sample]
    row_data = [sample] + [counts.get(bin_name, 0) for bin_name in bins]
    result_df.loc[row_idx] = row_data
    row_idx += 1

# Add total column
result_df['Total'] = result_df[bins].sum(axis=1)

print(f"\\nProcessed {len(sample_counts)} samples")

# Write the formatted table to file - two versions:
# 1. Tab-separated for Excel import
with open("combined_summary_minor_variants.tsv", 'w') as f:
    # Write header
    f.write("# Summary of Minor APOBEC3 Variants Across Samples\\n")
    f.write("# Each cell shows the number of variants in the given frequency range for each sample\\n")
    f.write("#\\n")
    
    # Write header row with tab separators
    f.write("Sample\\t" + "\\t".join(bins) + "\\tTotal\\n")
    
    # Write data rows with tab separators
    for _, row in result_df.iterrows():
        f.write(f"{row['Sample']}\\t")
        for bin_name in bins:
            f.write(f"{row[bin_name]}\\t")
        f.write(f"{row['Total']}\\n")

# 2. Nicely aligned text version for viewing in terminal
with open("combined_summary_minor_variants.txt", 'w') as f:
    # Write header
    f.write("# Summary of Minor APOBEC3 Variants Across Samples\\n")
    f.write("# Each cell shows the number of variants in the given frequency range for each sample\\n")
    f.write("#\\n")
    
    # Get column widths for proper alignment
    col_widths = {col: max(len(str(col)), max(len(str(val)) for val in result_df[col]) if len(result_df) > 0 else 0) + 2 
                  for col in result_df.columns}
    
    # Write header row with padding
    for col in result_df.columns:
        f.write(f"{col:<{col_widths[col]}}")
    f.write("\\n")
    
    # Write separator line
    total_width = sum(col_widths.values())
    f.write("-" * total_width + "\\n")
    
    # Write data rows with padding
    for _, row in result_df.iterrows():
        for col in result_df.columns:
            f.write(f"{row[col]:<{col_widths[col]}}")
        f.write("\\n")

print(f"Created combined summary with {len(result_df)} samples")
print(f"Summary files created in two formats:")
print(f"  - Tab-separated for Excel: combined_summary_minor_variants.tsv")
print(f"  - Formatted text for viewing: combined_summary_minor_variants.txt")
EOF

    # Make the script executable and run it
    chmod +x combine_summary.py
    python combine_summary.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
} 