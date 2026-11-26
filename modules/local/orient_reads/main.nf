process ORIENT_READS {
    // Reverse-complement single-end reads whose rRNA hits are all on the minus strand, so every read runs
    // 5'->3' along the gene. Long-read libraries mix both strands, and primer identification looks for primers
    // at read starts only.
    tag "$meta.id"
    label 'process_low'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.9.0--h9ee0642_0' :
        'biocontainers/seqkit:2.9.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads), path(cmsearch_tblout)

    output:
    tuple val(meta), path("${meta.id}.oriented.fastq.gz"), emit: reads
    path "versions.yml"                                  , emit: versions

    script:
    """
    # tblout column 10 is the strand. Reads with hits on both strands (chimeras) are left as they are.
    awk '!/^#/ { s[\$1] = s[\$1] \$10 } END { for (r in s) if (s[r] !~ /\\+/) print r }' ${cmsearch_tblout} > minus_ids.txt

    if [ -s minus_ids.txt ]; then
        seqkit grep -f minus_ids.txt ${reads} | seqkit seq -r -p -t dna > minus.fastq
        seqkit grep -v -f minus_ids.txt ${reads} > plus.fastq
        cat plus.fastq minus.fastq | gzip > ${meta.id}.oriented.fastq.gz
    else
        cp ${reads} ${meta.id}.oriented.fastq.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
