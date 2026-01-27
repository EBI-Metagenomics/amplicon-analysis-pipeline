
process PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS {
    tag "$meta.id"
    label 'light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(cmsearch_deoverlap_out), path(concat_primers_fasta)

    output:
    tuple val(meta), path("*primer_validation.tsv"), emit: primer_validation_out
    tuple val(meta), path("*primers.fasta")        , emit: validated_primers
    path "versions.yml"                            , emit: versions

    script:
    def se_flag = meta.single_end ? "--se" : ""
    """
    primer_val_classification -i $cmsearch_deoverlap_out -f $concat_primers_fasta -s ${meta.id} ${se_flag}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
