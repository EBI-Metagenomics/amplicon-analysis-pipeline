
process ITS_SANITY_CHECKER {
    tag "$meta.id"
    label 'light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(its_reads), path(itsonedb_assignments), path(unite_assignments)

    output:
    tuple val(meta), path("*.json"), emit: its_sanity_check_out
    tuple val(meta), path("*.tsv"), emit: its_sanity_check_out_tsv
    path "versions.yml"           , emit: versions

    script:
    """
    its_sanity_checker.py --itsonedb_output ${itsonedb_assignments} --unite_output ${unite_assignments} -r ${its_reads} -p ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
