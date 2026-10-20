
process CUTADAPT {
    tag "$meta.id"
    label 'very_light'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.6--py39hf95cd2a_1':
        'biocontainers/cutadapt:4.6--py39hf95cd2a_1' }"

    input:
    tuple val(meta), path(reads), path(primers)

    output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')          , emit: log
    tuple val(meta), path('*.json')         , emit: json
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // primers convention: [fwd, rev, fwd_rc, rev_rc], any of which may be empty
    def (fwd_primer, rev_primer, fwd_primer_rc, rev_primer_rc) = primers ?: []

    // Pass 1 trims each read's own primer from its 5' end, pass 2 the other
    // primer (reverse-complemented) from its 3' end, left by read-through
    def first_primer_args  = fwd_primer?.size() ? "-g file:${fwd_primer}" : ''
    def second_primer_args = rev_primer_rc?.size() ? "-a file:${rev_primer_rc}" : ''
    if (!meta.single_end) {
        first_primer_args  += rev_primer?.size() ? " -G file:${rev_primer}" : ''
        second_primer_args += fwd_primer_rc?.size() ? " -A file:${fwd_primer_rc}" : ''
    }

    def suffixes = meta.single_end ? [""] : ["_1", "_2"]
    def first_pass = suffixes.collect { s -> "${prefix}${s}.first_pass.fastq.gz" }
    def trimmed    = suffixes.collect { s -> "${prefix}${s}.trim.fastq.gz" }
    def first_pass_out = meta.single_end ? "-o ${first_pass[0]}" : "-o ${first_pass[0]} -p ${first_pass[1]}"
    def trimmed_out    = meta.single_end ? "-o ${trimmed[0]}" : "-o ${trimmed[0]} -p ${trimmed[1]}"

    if (!first_primer_args && !second_primer_args) {
        """
        touch ${prefix}.first_pass.cutadapt.log ${prefix}.second_pass.cutadapt.log
        touch ${prefix}.first_pass.cutadapt.json ${prefix}.second_pass.cutadapt.json
        touch ${trimmed.join(' ')}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
    else {
        """
        cutadapt \\
            --cores $task.cpus \\
            $args \\
            $first_pass_out \\
            $first_primer_args \\
            $reads \\
            --json ${prefix}.first_pass.cutadapt.json \\
            > ${prefix}.first_pass.cutadapt.log

        cutadapt \\
            --cores $task.cpus \\
            $args \\
            $trimmed_out \\
            $second_primer_args \\
            ${first_pass.join(' ')} \\
            --json ${prefix}.second_pass.cutadapt.json \\
            > ${prefix}.second_pass.cutadapt.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def trimmed = meta.single_end ? "${prefix}.trim.fastq.gz" : "${prefix}_1.trim.fastq.gz ${prefix}_2.trim.fastq.gz"
    """
    touch ${prefix}.first_pass.cutadapt.log ${prefix}.second_pass.cutadapt.log
    touch ${prefix}.first_pass.cutadapt.json ${prefix}.second_pass.cutadapt.json
    touch ${trimmed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
