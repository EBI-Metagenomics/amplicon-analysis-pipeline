process BBMAP_REFORMAT_STANDARDISE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(reads)
    val(interleaved)
    val(out_fmt)

    output:
    tuple val(meta), path("*_reformated.${out_fmt}")                                             , emit: reformated
    tuple val(meta), path("${prefix}_singleton.${out_fmt}")                                      , optional: true, emit: singleton
    tuple val("${task.process}"), val('bbmap'), eval('bbversion.sh | grep -v "Duplicate cpuset"'), emit: versions_bbmap, topic: versions
    path  "versions.yml"                                                                         , emit: versions
    path  "*.log"                                                                                , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def interleaved_args = task.ext.interleaved_args ?: ''
    def paired_args = task.ext.paired_args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def in_reads = (meta.single_end || interleaved) ? "in=${reads[0]}" : "in=${reads[0]} in2=${reads[1]}"
    def out_reads = meta.single_end ? "out=${prefix}_reformated.${out_fmt}" : "out=${prefix}_1_reformated.${out_fmt} out2=${prefix}_2_reformated.${out_fmt} outs=${prefix}_singleton.${out_fmt}"
    def interleaved_cmd = interleaved ? "int=t " + interleaved_args : ""
    def paired_cmd = meta.single_end ? "" : paired_args

    """
    maxmem=\$(echo \"$task.memory\"| sed 's/ GB/g/g')
    reformat.sh \\
        -Xmx\$maxmem \\
        $in_reads \\
        $out_reads \\
        $interleaved_cmd \\
        $paired_cmd \\
        threads=${task.cpus} \\
        ${args} \\
        &> ${prefix}.reformat.sh.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_1_reformated.${out_fmt}
    echo "" | gzip > ${prefix}_2_reformated.${out_fmt}
    echo "" | gzip > ${prefix}_singleton.${out_fmt}
    touch ${prefix}.reformat.sh.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
