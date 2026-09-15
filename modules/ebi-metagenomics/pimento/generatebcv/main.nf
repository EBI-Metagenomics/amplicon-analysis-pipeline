
process PIMENTO_GENERATEBCV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mi-pimento:1.0.2--pyhdfd78af_0':
        'biocontainers/mi-pimento:1.0.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(fwd_flag), val(rev_flag), path(fastq)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def var_region = "${meta.var_region}" ?: ""
    def assess_mcp_prop_prefix = "${prefix}_${var_region}"
    def strands = ""

    if (fwd_flag == "auto" && rev_flag == "auto") {
        strands = "FR"
    } else if (fwd_flag == "auto") {
        strands = "F"
    } else if (rev_flag == "auto") {
        strands = "R"
    }

    if (strands == "") {
        """
        touch ${assess_mcp_prop_prefix}_mcp_cons.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mi-pimento: \$( pimento --version | cut -d" " -f3 )
        END_VERSIONS
        """
    } else {
        """
        pimento \\
            gen_bcv \\
            -i ${fastq} \\
            -st ${strands} \\
            -o ${assess_mcp_prop_prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mi-pimento: \$( pimento --version | cut -d" " -f3 )
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def var_region = "${meta.var_region}"
    def assess_mcp_prop_prefix = "${prefix}_${var_region}"

    """
    touch ${assess_mcp_prop_prefix}_mcp_cons.tsv

    echo 'dummy' > ${assess_mcp_prop_prefix}_mcp_cons.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mi-pimento: \$( pimento --version | cut -d" " -f3 )
    END_VERSIONS
    """
}
