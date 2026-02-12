
process MAPSEQ2BIOM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:1.4.12--pyhdfd78af_0' :
        'biocontainers/mgnify-pipelines-toolkit:1.4.12--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(msq), path(db_otu), val(db_label)

    output:
    tuple val(meta), path("${meta.id}_${task.ext.db_label}.txt"), emit: krona_input
    tuple val(meta), path("${meta.id}_${task.ext.db_label}.tsv"), emit: biom_out
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mapseq2biom \
        ${args} \
        --krona ${prefix}_${task.ext.db_label}.txt \
        --no-tax-id-file ${prefix}.notaxid.tsv \
        --label ${db_label} \
        --query ${msq} \
        --otu-table ${db_otu} \
        --out-file ${prefix}_${task.ext.db_label}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: \$(get_mpt_version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_${task.ext.db_label}.txt
    touch ${prefix}_${task.ext.db_label}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapseq2biom: 0.1.1
    END_VERSIONS
    """
}
