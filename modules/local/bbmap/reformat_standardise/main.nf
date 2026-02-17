process BBMAP_REFORMAT_STANDARDISE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5a/5aae5977ff9de3e01ff962dc495bfa23f4304c676446b5fdf2de5c7edfa2dc4e/data' :
        'community.wave.seqera.io/library/bbmap_pigz:07416fe99b090fa9' }"

    input:
    tuple val(meta), path(reads)
    val(out_fmt)

    output:
    tuple val(meta), path("*_reformated.${out_fmt}")       , emit: reformated
    tuple val(meta), path("${prefix}_singleton.${out_fmt}"), optional: true, emit: singleton
    path  "versions.yml"                                      , emit: versions
    path  "*.log"                                             , emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    single_file = (reads instanceof Collection) ? (reads.size() == 1) : true
    in_reads  = single_file ? "in=${reads[0]}" : "in=${reads[0]} in2=${reads[1]}"
    out_reads = meta.single_end ? "out=${prefix}_reformated.${out_fmt}" : "out=${prefix}_1_reformated.${out_fmt} out2=${prefix}_2_reformated.${out_fmt} outs=${prefix}_singleton.${out_fmt}"
    interleaved_cmd = meta.interleaved ? "int=t verifyinterleaved=t" : ""
    paired_cmd = meta.single_end ? "" : "addslash=t spaceslash=t verifypaired=t"

    """
    maxmem="${task.memory.toGiga()}G"
    reformat.sh \\
        -Xmx\$maxmem \\
        $in_reads \\
        $out_reads \\
        $interleaved_cmd \\
        $paired_cmd \\
        threads=${task.cpus} \\
        allowidenticalnames=t \\
        trimreaddescription=t \\
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
    touch ${prefix}.repair.sh.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh | grep -v "Duplicate cpuset")
    END_VERSIONS
    """
}
