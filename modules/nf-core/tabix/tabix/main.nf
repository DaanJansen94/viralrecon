process TABIX_TABIX {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(tab)

    output:
    tuple val(meta), path("*.tbi"), emit: tbi
    tuple val(meta), path(tab), emit: gz
    tuple val(meta), path("*.csi"), optional:true, emit: csi
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    
    """
    # First check if the input file exists and has content
    if [ ! -s "${tab}" ]; then
        echo "Error: Input file ${tab} is empty or does not exist"
        exit 1
    fi
    
    # Debug information
    echo "Input file: ${tab}"
    ls -l ${tab}
    
    # Check if input is already compressed
    if [[ "${tab}" != *.gz ]]; then
        echo "Input file is not compressed"
        # Compress the input file with bgzip
        bgzip -c ${tab} > ${tab}.gz
        # Index the compressed file
        tabix $args ${tab}.gz
    else
        echo "Input file is already compressed"
        # Index the compressed file
        tabix $args ${tab}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${tab}.tbi
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
