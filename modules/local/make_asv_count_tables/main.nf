
process MAKE_ASV_COUNT_TABLES {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(maps), path(reads), path(filter_list)

    output:
    tuple val(meta), path("*_asv_read_counts.tsv"), optional: true, emit: asv_read_counts
    path "versions.yml"                           , emit: versions

    script:
    """
    # Passing an empty taxonomy table to make_asv_count_table means that it does not filter by taxonomy
    printf 'ASV\\tSuperkingdom\\tKingdom\\tPhylum\\tClass\\tOrder\\tFamily\\tGenus\\tSpecies\\n' > empty_tax.tsv

    if [[ ${meta.single_end} = true ]]; then
        zcat ${reads} | awk 'NR % 4 == 1' > headers.txt
        make_asv_count_table -t empty_tax.tsv -f ${maps} -a ${filter_list} -hd headers.txt -s ${meta.id}_${meta.var_region}
    else
        zcat ${reads[0]} | awk 'NR % 4 == 1' > headers.txt
        make_asv_count_table -t empty_tax.tsv -f ${maps[0]} -a ${filter_list} -hd headers.txt -s ${meta.id}_${meta.var_region}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${meta.var_region}_asv_read_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
