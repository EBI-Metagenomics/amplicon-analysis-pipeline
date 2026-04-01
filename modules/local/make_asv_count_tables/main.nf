
process MAKE_ASV_COUNT_TABLES {
    tag "$meta.id"
    label 'medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(maps), path(asvtaxtable), path(reads), path(extracted_var_path), val(db_label)

    output:
    tuple val(meta), path("*asv_krona_counts.txt"), optional: true, emit: asv_krona_counts
    tuple val(meta), path("*asv_read_counts.tsv") , optional: true, emit: asv_read_counts
    path "versions.yml"                           , emit: versions

    script:
    """
    if [[ ${meta.single_end} = true ]]; then
        zcat $reads | awk 'NR % 4 == 1' > headers.txt
        make_asv_count_table -t $asvtaxtable -f $maps -a $extracted_var_path -hd ./headers.txt  -s ${meta.id}_${meta.var_region}_${db_label}
    else
        zcat ${reads[0]} | awk 'NR % 4 == 1' > headers.txt
        make_asv_count_table -t $asvtaxtable -f ${maps[0]} -a $extracted_var_path -hd ./headers.txt  -s ${meta.id}_${meta.var_region}_${db_label}
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
