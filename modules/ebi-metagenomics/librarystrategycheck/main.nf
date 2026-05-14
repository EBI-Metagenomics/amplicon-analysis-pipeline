
process LIBRARYSTRATEGYCHECK {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.2.1--pyhdfd78af_0':
        'biocontainers/mgnify-pipelines-toolkit:1.2.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bcv)

    output:
    tuple val(meta), env('check_out'), emit: library_check_out
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    library_strategy_check \\
        -i ${bcv} \\
        -s ${prefix} \\
        -o ./

    check_out=\$(cat ${prefix}_library_check_out.txt)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    """
    check_out="AMPLICON"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
