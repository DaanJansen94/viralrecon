process APOBEC3_ANALYSIS_IVAR {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9.5 conda-forge::pandas=1.3.5 bioconda::pysam=0.19.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9.5' :
        'python:3.9.5' }"

    publishDir "${params.outdir}/variants/ivar/apobec3", mode: params.publish_dir_mode, pattern: "*_apobec_mutations.csv"
    publishDir "${params.outdir}/variants/ivar/apobec3", mode: params.publish_dir_mode, pattern: "apobec_summary.txt"
    publishDir "${params.outdir}/variants/ivar/apobec3", mode: params.publish_dir_mode, pattern: "read_counts.csv"
    publishDir "${params.outdir}/variants/ivar/apobec3", mode: params.publish_dir_mode, pattern: "versions.yml"

    input:
    tuple val(meta), path(variants)
    path(reference)
    path(apobec3_script)
    path(count_reads_script)
    path(samplesheet)

    output:
    tuple val(meta), path("*_apobec_mutations.csv"), emit: csv
    path("apobec_summary.txt"), emit: summary
    path("read_counts.csv"), emit: read_counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p /tmp/package_install
    export PIP_TARGET=/tmp/package_install
    export PYTHONPATH="/tmp/package_install"
    pip install --no-cache-dir --target=/tmp/package_install numpy==1.22.4 pandas==1.3.5 pysam==0.19.0

    # Count raw reads from samplesheet
    python ${count_reads_script} ${samplesheet}

    # Run APOBEC3 analysis
    python ${apobec3_script} ${variants} ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}