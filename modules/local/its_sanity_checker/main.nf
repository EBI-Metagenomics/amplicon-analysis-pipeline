
process ITS_SANITY_CHECKER {
    tag "$meta.id"
    label 'light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), val(read_assignments)

    output:
    tuple val(meta), path("*_its_sanity_check.json")   , emit: its_sanity_check_out
    tuple val(meta), path("*_its_sanity_check_mqc.tsv"), emit: its_sanity_check_out_mqc
    tuple val(meta), path("read_assignments.json")     , emit: read_assignment_counts
    path "versions.yml"                                , emit: versions

    script:
    def serializable = read_assignments.collectEntries { k, v -> [k, v.toString()] }
    def read_assignments_json = new groovy.json.JsonBuilder(serializable).toString()
    """
    echo '${read_assignments_json}' > read_assignments.json
    its_sanity_checker.py --read_assignments read_assignments.json -p ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
