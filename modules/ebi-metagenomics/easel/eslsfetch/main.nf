
process EASEL_ESLSFETCH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::easel=0.49"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/easel:0.49--h031d066_0':
        'biocontainers/easel:0.49--h031d066_0' }"

    input:
    tuple val(meta), path(fasta), path(cmsearch_deoverlap_out)

    output:
    tuple val(meta), path("*easel_coords.fasta"),           emit: easel_coords
    tuple val(meta), path("*matched_seqs_with_coords.txt"), emit: matched_seqs_with_coords
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def is_compressed = fasta.getExtension() == "gz"
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta

    """

    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    awk \\
        '{print \$1"-"\$3"/"\$8"-"\$9" "\$8" "\$9" "\$1}' \\
        $cmsearch_deoverlap_out \\
        > ${prefix}.matched_seqs_with_coords.txt

    esl-sfetch \\
        --index \\
        $fasta_name

    esl-sfetch \\
        -Cf \\
        $fasta_name \\
        ${prefix}.matched_seqs_with_coords.txt \\
        > ${prefix}_easel_coords.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        easel: \$(esl-sfetch -h | grep -o '^# Easel [0-9.]*' | sed 's/^# Easel *//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.matched_seqs_with_coords.txt
    touch ${prefix}_easel_extracted.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        easel: \$(esl-sfetch -h | grep -o '^# Easel [0-9.]*' | sed 's/^# Easel *//')
    END_VERSIONS
    """
}
