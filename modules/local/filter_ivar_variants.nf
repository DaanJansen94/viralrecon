process FILTER_IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::pandas=1.5.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9.5' :
        'python:3.9.5' }"

    publishDir "${params.outdir}/variants/ivar", mode: params.publish_dir_mode, pattern: "*.filtered.tsv"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.filtered.tsv"), emit: tsv
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
    df = pd.read_csv("${tsv}", sep='\t')
    
    # Store original counts
    original_count = len(df)
    
    # Apply filtering based on ALT_DP
    filtered_df = df[df['ALT_DP'] >= ${min_alt_dp}]
    
    # Optionally filter out indels (keep only SNPs)
    if ${filter_indels}:
        # Keep only SNPs (where REF and ALT are single bases)
        filtered_df = filtered_df[
            (filtered_df['REF'].str.len() == 1) & 
            (filtered_df['ALT'].str.len() == 1) &
            (filtered_df['REF'].str.match("^[ATGC]\$")) &
            (filtered_df['ALT'].str.match("^[ATGC]\$"))
        ]
    
    # Write the filtered data to a new file
    filtered_df.to_csv("${prefix}.filtered.tsv", sep='\t', index=False)
    
    print(f"Processed {original_count} variants")
    print(f"After filtering with ALT_DP >= ${min_alt_dp} and indel removal: ${filter_indels}")
    print(f"Remaining variants: {len(filtered_df)}")
    
except Exception as e:
    print(f"Error processing file: {str(e)}", file=sys.stderr)
    sys.exit(1)
EOF

    cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
    pandas: \$(python -c "import pandas; print(pandas.__version__)")
END_VERSIONS
    """
} 