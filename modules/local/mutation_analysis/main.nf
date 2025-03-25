process MUTATION_ANALYSIS_IVAR {
    tag "all_samples"
    label 'process_low'

    conda "conda-forge::python=3.9.5 conda-forge::pandas=1.3.5 bioconda::pysam=0.19.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9.5' :
        'python:3.9.5' }"

    publishDir "${params.outdir}/variants/ivar/mutation_analysis", mode: params.publish_dir_mode, pattern: "mutations_and_apobec3_summary.txt"
    publishDir "${params.outdir}/variants/ivar/mutation_analysis", mode: params.publish_dir_mode, pattern: "mutations_and_apobec3_summary.csv"
    publishDir "${params.outdir}/variants/ivar/mutation_analysis", mode: params.publish_dir_mode, pattern: "versions.yml"

    input:
    tuple val(meta), path(variants)
    path(stats)
    path(mutation_script)

    output:
    path("mutations_and_apobec3_summary.txt"), emit: summary
    path("mutations_and_apobec3_summary.csv"), emit: csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p /tmp/package_install
    export PIP_TARGET=/tmp/package_install
    export PYTHONPATH="/tmp/package_install"
    pip install --no-cache-dir --target=/tmp/package_install numpy==1.22.4 pandas==1.3.5 pysam==0.19.0
    
    # Create directories for input files
    mkdir -p variants stats

    # Copy input files with their original names
    for variant_file in ${variants}; do
        cp -L "\$variant_file" variants/\$(basename "\$variant_file")
    done

    for stats_file in ${stats}; do
        cp -L "\$stats_file" stats/\$(basename "\$stats_file")
    done

    # Run mutation analysis script - now using local directories
    python ${mutation_script} variants stats

    cat <<-END_VERSIONS > versions.yml
    "NFCORE_VIRALRECON:ILLUMINA:MUTATION_ANALYSIS_IVAR":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
} 