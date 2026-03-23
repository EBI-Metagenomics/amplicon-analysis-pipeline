
process MAPSEQ {
    tag "$meta.id"
    label 'medium'

    conda "bioconda::mapseq=2.1.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mapseq:2.1.1b--h3ab3c3b_0':
        'biocontainers/mapseq:2.1.1b--h3ab3c3b_0' }"

    input:
    tuple val(meta), path(subunit_reads)
    tuple path(db_fasta), path(db_tax), path(db_mscluster)

    output:
    tuple val(meta), path("*.mseq"), emit: mseq
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mapseq \\
        -seed 12 \\
        $subunit_reads \\
        $db_fasta \\
        $db_tax \\
        -nthreads $task.cpus \\
        $args \\
        > ${prefix}.mseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapseq: \$(mapseq -h 2>&1 | head -1 | cut -d" " -f2 | cut -d"v" -f2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.mseq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mapseq: \$(mapseq -h 2>&1 | head -1 | cut -d" " -f2 | cut -d"v" -f2)
    END_VERSIONS
    """
}
