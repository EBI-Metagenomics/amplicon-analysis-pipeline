
process FASTQSUFFIXHEADERCHECK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.2.11--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:1.2.11--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.json"), emit: json
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq_input = meta.single_end ? "-f ${fastq}" : "-f ${fastq[0]} -r ${fastq[1]}"

    """
    fastq_suffix_header_check \\
    ${fastq_input} \\
    -s ${prefix} \\
    -o ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_suffix_header_err.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
