process FILTER_MAPSEQ_BY_MASK {
    tag "${meta.id}:${meta.db_label}:${meta.target}"
    label 'process_single'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(mseq), path(masked_reads)

    output:
    tuple val(meta), path("${meta.id}_${meta.db_label}_${meta.target}_filtered.mseq"), emit: mseq

    script:
    def filtered_mseq = "${meta.id}_${meta.db_label}_${meta.target}_filtered.mseq"
    """
    awk -v masked_fasta='${masked_reads}' '
        BEGIN {
            while ((getline line < masked_fasta) > 0) {
                if (line ~ /^>/) {
                    sub(/^>/, "", line)
                    sub(/[[:space:]].*\$/, "", line)
                    sub(/\\r\$/, "", line)
                    allowed[line] = 1
                }
            }
            close(masked_fasta)
        }

        /^#/ {
            print
            next
        }

        {
            query = \$1
            slash = index(query, "/")
            read_id = query

            if (slash) {
                prefix = substr(query, 1, slash - 1)

                # MAPseq appends -<subunit>/<start>-<end> to the original read
                # ID. Taking the last dash before the coordinates also supports
                # read IDs that themselves contain dashes.
                last_dash = 0
                for (i = 1; i <= length(prefix); i++) {
                    if (substr(prefix, i, 1) == "-")
                        last_dash = i
                }
                read_id = last_dash ? substr(prefix, 1, last_dash - 1) : prefix
            }

            if (allowed[read_id])
                print
        }
    ' '${mseq}' > '${filtered_mseq}'
    """
}

