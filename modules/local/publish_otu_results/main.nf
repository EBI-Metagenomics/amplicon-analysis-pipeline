
process PUBLISH_OTU_RESULTS {
    tag "$meta.id"
    label 'light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(mseq), path(krona_input), path(biom), path(html)

    output:
    tuple val(meta), path("*.mseq", includeInputs: true), emit: mseq
    tuple val(meta), path("*.txt", includeInputs: true) , emit: krona_input
    tuple val(meta), path("*.tsv", includeInputs: true) , emit: biom
    tuple val(meta), path("*.html", includeInputs: true), emit: html
    path "versions.yml"                                 , emit: versions

    script:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
