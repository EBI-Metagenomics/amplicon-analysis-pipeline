process MAPSEQ2BIOM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.4.12--pyhdfd78af_0' :
        'biocontainers/mgnify-pipelines-toolkit:1.4.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(msq)
    tuple path(db_otu), val(db_label)

    output:
    tuple val(meta), path("${prefix}.txt")        , emit: krona_input
    tuple val(meta), path("${prefix}.tsv")        , emit: biom_out
    tuple val(meta), path("${prefix}.notaxid.tsv"), emit: biom_notaxid_out, optional: true
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mapseq2biom \
        ${args} \
        --krona ${prefix}.txt \
        --no-tax-id-file ${prefix}.notaxid.tsv \
        --label ${db_label} \
        --query ${msq} \
        --otu-table ${db_otu} \
        --out-file ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.txt
    touch ${prefix}.notaxid.tsv
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """
}
