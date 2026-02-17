
process FILTER_MASKED_N {
    // Remove reads that have at least 10% N bases with seqkit commands
    tag "$meta.id"
    label 'process_medium'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_ITS_rRNA.fa"), emit: filtered_its_fasta, optional: true
    path "versions.yml"                   , emit: versions

    script:
    """ 
    seqkit fx2tab -B ACGT ${reads} \
        | awk '\$3 > 50' > temp_fasta.tab

    num_lines=\$(wc -l temp_fasta.tab | cut -d' ' -f1)

    if [ \$num_lines -gt 0 ]
    then
        seqkit tab2fx temp_fasta.tab -o ${meta.id}_ITS_rRNA.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit | sed '3!d; s/Version: //' )
    END_VERSIONS
    """
}
