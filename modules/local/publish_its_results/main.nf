
process PUBLISH_ITS_RESULTS {
    tag "$meta.id"
    label 'light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(results_files)

    output:
    tuple val(meta), path("*.mseq", includeInputs: true), emit: mapseq_output
    tuple val(meta), path("*.txt", includeInputs: true) , emit: krona_input
    tuple val(meta), path("*.tsv", includeInputs: true) , emit: biom_output
    tuple val(meta), path("*.html", includeInputs: true), emit: krona_output
    tuple val(meta), path("*.fa", includeInputs: true)  , emit: its_reads
    path "versions.yml"                                , emit: versions

    script:
    """
    echo "Publishing ITS results!"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
