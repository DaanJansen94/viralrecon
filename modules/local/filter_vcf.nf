process FILTER_VCF {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::tabix=0.2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:0.2.6--ha92aebf_0' :
        'quay.io/biocontainers/tabix:0.2.6--ha92aebf_0' }"

    publishDir "${params.outdir}/variants/ivar", mode: params.publish_dir_mode, pattern: "*.filter.vcf.gz"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.filter.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_alt_dp = params.min_alt_dp ?: 5
    def filter_indels = params.filter_indels ? "True" : "False"
    
    """
    # Extract the headers directly to the output file
    zcat ${vcf} | grep "^#" > ${prefix}.filter.vcf

    # Count the total number of variants for statistics
    TOTAL=\$(zcat ${vcf} | grep -v "^#" | wc -l)
    
    # Filter the variants based on ALT_DP and optionally remove indels
    if [ "${filter_indels}" = "True" ]; then
        zcat ${vcf} | grep -v "^#" | awk -F'\\t' '{
            # Extract reference and alternate alleles
            ref = \$4
            alt = \$5
            
            # Check if this is a SNP (same length ref and alt)
            is_snp = (length(ref) == length(alt))
            
            # Parse FORMAT field to find ALT_DP index
            split(\$9, fmt, ":")
            split(\$10, vals, ":")
            
            alt_dp = 0
            for (i=1; i<=length(fmt); i++) {
                if (fmt[i] == "ALT_DP") {
                    alt_dp = vals[i]
                    break
                }
            }
            
            # Only print variants that meet criteria
            if (alt_dp >= ${min_alt_dp} && is_snp) {
                print \$0
            }
        }' >> ${prefix}.filter.vcf
    else
        zcat ${vcf} | grep -v "^#" | awk -F'\\t' '{
            # Parse FORMAT field to find ALT_DP index
            split(\$9, fmt, ":")
            split(\$10, vals, ":")
            
            alt_dp = 0
            for (i=1; i<=length(fmt); i++) {
                if (fmt[i] == "ALT_DP") {
                    alt_dp = vals[i]
                    break
                }
            }
            
            # Only print variants that meet criteria
            if (alt_dp >= ${min_alt_dp}) {
                print \$0
            }
        }' >> ${prefix}.filter.vcf
    fi
    
    # Count how many variants were filtered
    REMAINING=\$(grep -v "^#" ${prefix}.filter.vcf | wc -l)
    FILTERED=\$((\$TOTAL - \$REMAINING))
    
    # Print statistics
    echo "Total variants: \$TOTAL"
    echo "Filtered variants: \$FILTERED"
    echo "Remaining variants: \$REMAINING"
    
    # Compress the filtered VCF with bgzip for tabix compatibility
    bgzip -f ${prefix}.filter.vcf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: \$(bgzip -h 2>&1 | head -n 1 | sed 's/^.*(htslib) //; s/ .*\$//')
        tabix: \$(tabix -h 2>&1 | head -n 1 | sed 's/^.*(htslib) //; s/ .*\$//')
    END_VERSIONS
    """
} 