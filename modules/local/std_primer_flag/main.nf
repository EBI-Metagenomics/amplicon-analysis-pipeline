
process STD_PRIMER_FLAG {
    // Check for presence of standard library of primers (stored in ./data/standard_primers)
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mi-pimento:${params.pimento_version}":
        "biocontainers/mi-pimento:${params.pimento_version}"}"
    input:
    tuple val(meta), path(reads_merged)
    path(std_primer_library)

    output:
    tuple val(meta), path("*std_primers.fasta"), emit: std_primer_out
    path "*std_primer_out.txt"
    path "versions.yml"                        , emit: versions

    script:
    def std_primer_library_arg = "${std_primer_library}" ? "-p ${std_primer_library}" : ""

    """
    pimento std -i ${reads_merged} ${std_primer_library_arg} --merged -o ${meta.id}_${meta.var_region}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mi-pimento: \$( pimento --version | cut -d" " -f3 )
    END_VERSIONS
    """

}
