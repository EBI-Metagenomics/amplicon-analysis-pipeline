
process CLASSIFY_VAR_REGIONS {
    tag "$meta.id"
    label 'light'
    conda "${moduleDir}/../../../conf/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(cmsearch_deoverlap_out)

    output:
    tuple val(meta), path ("*S.V*.txt", arity: '0..2')          , optional: true, emit: classify_var_regions
    tuple val(meta), val("concat"), path("*.concat.regions.txt"), optional: true, emit: concat_var_regions
    tuple val(meta), path("*.tsv")                              , optional: true, emit: cdch_out
    path "versions.yml"                                         , emit: versions

    script:
    """
    classify_var_regions -d ./ -o ${meta.id} --statistics $cmsearch_deoverlap_out
    
    num_files="\$(find . -name '*S.V*.txt' | wc -l  )"
    if [ \$num_files -gt 1 ]; then
        cat *S.V*.txt > ${meta.id}.concat.regions.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

}
