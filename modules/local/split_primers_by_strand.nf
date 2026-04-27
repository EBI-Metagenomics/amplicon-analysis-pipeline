
process SPLIT_PRIMERS_BY_STRAND {
    // Split the concatenated primer files into two, one containing the fwd primers, and one for the rev
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"
    
    input:
    tuple val(meta), path(concat_primers)

    output:
    tuple val(meta), path("fwd_primer.fasta"), path("rev_primer.fasta"), emit: stranded_primer_out
    path "versions.yml"                                                , emit: versions

    script:
    """
    if [[ -s ${concat_primers} ]]; then

        fwd_flag="\$(grep -A 1 '^>.*F\$' ${concat_primers} | grep -v -- '^--\$' || true)"
        rev_flag="\$(grep -A 1 '^>.*R\$' ${concat_primers} | grep -v -- '^--\$' || true)"

        if [[ ! -z \$fwd_flag ]]; then
            grep -A 1 '^>.*F\$' ${concat_primers} | grep -v -- '^--\$' > ./fwd_primer.fasta
        else
            touch ./fwd_primer.fasta
        fi

        if [[ ! -z \$rev_flag ]]; then
            grep -A 1 '^>.*R\$' ${concat_primers} | grep -v -- '^--\$' > ./rev_primer.fasta
        else
            touch ./rev_primer.fasta
        fi
    else
        touch ./fwd_primer.fasta
        touch ./rev_primer.fasta
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
