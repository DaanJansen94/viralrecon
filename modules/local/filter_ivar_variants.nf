process FILTER_IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pandas=1.5.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9.5' :
        'python:3.9.5' }"

    publishDir "${params.outdir}/variants/ivar", mode: params.publish_dir_mode, pattern: "*.filtered.tsv"
    publishDir "${params.outdir}/variants/ivar/log", mode: params.publish_dir_mode, pattern: "*.filter_stats.log"
    publishDir "${params.outdir}/variants/ivar/masked", mode: params.publish_dir_mode, pattern: "*.masked.tsv"

    input:
    tuple val(meta), path(tsv)
    path mask_file

    output:
    tuple val(meta), path("*.filtered.tsv"), emit: tsv
    path "*.filter_stats.log", emit: log
    path "*.masked.tsv", emit: masked
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_alt_dp = params.min_alt_dp ?: 5
    def filter_indels = params.filter_indels ? "True" : "False"
    """
    mkdir -p /tmp/package_install
    export PIP_TARGET=/tmp/package_install
    export PYTHONPATH="/tmp/package_install"
    pip install --no-cache-dir --target=/tmp/package_install numpy==1.22.4 pandas==1.3.5

    python3 << EOF
import sys
import pandas as pd

# Read the TSV file
try:
    df = pd.read_csv("${tsv}", sep='\\t')
    
    # Store original counts
    original_count = len(df)
    print(f"Total variants before filtering: {original_count}")
    
    # Track reasons for filtering
    filter_stats = {
        'alt_dp': 0,
        'indels': 0,
        'masked': 0
    }
    
    # Apply filtering based on ALT_DP
    alt_dp_filtered = df[df['ALT_DP'] < ${min_alt_dp}]
    filter_stats['alt_dp'] = len(alt_dp_filtered)
    df = df[df['ALT_DP'] >= ${min_alt_dp}]
    print(f"Variants filtered by ALT_DP < ${min_alt_dp}: {filter_stats['alt_dp']}")
    
    # Optionally filter out indels (keep only SNPs)
    if ${filter_indels}:
        # Keep only SNPs (where REF and ALT are single bases)
        indel_filtered = df[
            ~((df['REF'].str.len() == 1) & 
            (df['ALT'].str.len() == 1) &
            (df['REF'].str.match("^[ATGC]\$")) &
            (df['ALT'].str.match("^[ATGC]\$")))
        ]
        filter_stats['indels'] = len(indel_filtered)
        df = df[
            (df['REF'].str.len() == 1) & 
            (df['ALT'].str.len() == 1) &
            (df['REF'].str.match("^[ATGC]\$")) &
            (df['ALT'].str.match("^[ATGC]\$"))
        ]
        print(f"Variants filtered as indels: {filter_stats['indels']}")
    
    # Apply masking if mask file is provided
    if "${mask_file}" != "None":
        mask_df = pd.read_csv("${mask_file}")
        # Create a list of positions to mask
        masked_positions = []
        for _, row in mask_df.iterrows():
            masked_positions.extend(range(row['Minimum'], row['Maximum'] + 1))
        
        # Filter out variants in masked positions
        masked_filtered = df[df['POS'].isin(masked_positions)]
        filter_stats['masked'] = len(masked_filtered)
        
        # Save masked variants to a separate file with additional information
        if len(masked_filtered) > 0:
            masked_filtered['MASKED_REGION'] = masked_filtered['POS'].apply(
                lambda x: next((f"{row['Minimum']}-{row['Maximum']}" 
                              for _, row in mask_df.iterrows() 
                              if row['Minimum'] <= x <= row['Maximum']), '')
            )
            masked_filtered.to_csv("${prefix}.masked.tsv", sep='\\t', index=False)
        
        df = df[~df['POS'].isin(masked_positions)]
        print(f"Variants filtered by masking: {filter_stats['masked']}")
    
    # Write the filtered data to a new file
    df.to_csv("${prefix}.filtered.tsv", sep='\\t', index=False)
    
    # Write detailed stats to log file
    with open("${prefix}.filter_stats.log", "w") as f:
        f.write("Variant Filtering Statistics\\n")
        f.write("==========================\\n")
        f.write(f"Total variants before filtering: {original_count}\\n")
        f.write(f"Variants filtered by ALT_DP < ${min_alt_dp}: {filter_stats['alt_dp']}\\n")
        if ${filter_indels}:
            f.write(f"Variants filtered as indels: {filter_stats['indels']}\\n")
        if "${mask_file}" != "None":
            f.write(f"Variants filtered by masking: {filter_stats['masked']}\\n")
        f.write(f"Remaining variants: {len(df)}\\n")
    
    print(f"Processed {original_count} variants")
    print(f"After filtering with ALT_DP >= ${min_alt_dp} and indel removal: ${filter_indels}")
    print(f"Remaining variants: {len(df)}")
    
except Exception as e:
    print(f"Error processing file: {str(e)}", file=sys.stderr)
    sys.exit(1)
EOF

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    pandas: \$(python -c "import pandas; print(pandas.__version)")
END_VERSIONS
    """
} 