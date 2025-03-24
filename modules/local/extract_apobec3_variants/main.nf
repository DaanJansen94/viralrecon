process EXTRACT_APOBEC3_VARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.9.5 conda-forge::pandas=1.3.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9.5' :
        'python:3.9.5' }"

    publishDir "${params.outdir}/variants/ivar/apobec3", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(variants)
    path(extract_script)
    val(min_freq)
    val(min_depth)

    output:
    tuple val(meta), path("*.apobec3.tsv"), emit: tsv
    tuple val(meta), path("*.apobec3.major.tsv"), emit: major_tsv
    tuple val(meta), path("*.apobec3.minor.tsv"), emit: minor_tsv
    path "combined_summary_minor_variants.tsv", optional: true, emit: summary
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Set up temporary directory for pip packages
    mkdir -p /tmp/package_install
    export PIP_TARGET=/tmp/package_install
    export PYTHONPATH="/tmp/package_install"
    
    # Install numpy first with specific version compatible with pandas 1.3.5
    pip install --no-cache-dir --target=/tmp/package_install numpy==1.21.6
    
    # Then install pandas
    pip install --no-cache-dir --target=/tmp/package_install pandas==1.3.5
    
    mkdir -p variants out

    # Copy input files to variants directory
    cp ${variants} variants/

    # Run the extraction script
    python ${extract_script} \\
        --input_dir variants/ \\
        --output_dir out/ \\
        --min_freq ${min_freq} \\
        --min_depth ${min_depth}

    # Copy output to current directory
    cp out/*.apobec3*.tsv .
    if [ -f out/combined_summary_minor_variants.tsv ]; then
        cp out/combined_summary_minor_variants.tsv .
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
} 