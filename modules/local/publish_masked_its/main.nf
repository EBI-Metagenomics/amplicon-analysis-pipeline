process PUBLISH_MASKED_ITS {
    tag "$meta.id"

    input:
    tuple val(meta), path(masked_fasta)

    output:
    tuple val(meta), path("*_ITS_rRNA.fa", includeInputs: true), emit: masked_out

    script:
    """
    true
    """
}
