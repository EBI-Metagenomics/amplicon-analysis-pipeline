
process EXTRACT_ASVS_LEFT {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"
    
    input:
    tuple val(meta), path(asvs_left), path(fasta), path(asvtaxtable), val(db_label)

    output:
    tuple val(meta), path("*asv_seqs.fasta"), emit: asv_seqs_out
    tuple val(meta), path("*asv_tax.tsv")   , emit: asv_taxa_out
    path "versions.yml"                     , emit: versions

    script:
    """
    cat ${asvs_left} | cut -f1 > concatenated_asvs_left.txt

    grep -A 1 -wFf concatenated_asvs_left.txt ${fasta} > ${meta.id}_asv_seqs_temp.fasta
    grep -v '-' ${meta.id}_asv_seqs_temp.fasta > ${meta.id}_asv_seqs.fasta

    head -1 ${asvtaxtable} > ${meta.id}_${db_label}_asv_tax.tsv # keep the header
    grep -wFf concatenated_asvs_left.txt ${asvtaxtable} >> ${meta.id}_${db_label}_asv_tax.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
