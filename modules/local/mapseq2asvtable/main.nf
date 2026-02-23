
process MAPSEQ2ASVTABLE {
    tag "$meta.id"
    label 'very_light' // Will likely need to give this task more CPUs 
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(mapseq_out), val(db_label)

    output:
    tuple val(meta), path("*.tsv"), emit: asvtaxtable
    path "versions.yml"           , emit: versions

    script:
    """
    mapseq_to_asv_table -i $mapseq_out -l $db_label -s ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
