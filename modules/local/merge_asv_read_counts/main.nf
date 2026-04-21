
process MERGE_ASV_READ_COUNTS {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(mapseq_asv_counts, stageAs: "?/*")

    output:
    tuple val(meta), path("*asv_read_counts.tsv"), emit: merged_counts
    path "versions.yml"                          , emit: versions

    script:
    """
    cat ${mapseq_asv_counts} | sort | uniq > ${meta.id}_${meta.var_region}_asv_read_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
