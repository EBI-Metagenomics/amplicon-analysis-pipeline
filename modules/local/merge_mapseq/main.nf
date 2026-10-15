process MERGE_MAPSEQ {
    tag "${meta.id}:${meta.db_label}"
    label 'light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    // Use a unique staged pattern because the input MAPseq files have the
    // same original filename. Target-labelled symlinks are created below.
    // For a collection containing only one file, the staged filename becomes mseq_1.mseq
    // For multiple files, the names become mseq_1.mseq, mseq_2.mseq, ...
    tuple val(meta), path(mseq_files, stageAs: 'mseq_?.mseq')

    output:
    tuple val(meta), path("${meta.id}_${meta.db_label}.mseq"), emit: mseq

    script:
    def merged_mseq = "${meta.id}_${meta.db_label}.mseq"
    def staged_mseq_files = mseq_files instanceof List ? mseq_files : [mseq_files]
    def targets = meta.targets instanceof Collection ? meta.targets : [meta.targets]
    def target_symlinks = staged_mseq_files.withIndex().collect { mseq, index ->
        def target = targets[index] ?: "target_${index + 1}"
        def symlink = "${meta.id}_${meta.db_label}_${target}_${index + 1}.mseq"
        "ln -s -- '${mseq}' '${symlink}'"
    }.join('\n')
    def input_files = staged_mseq_files.withIndex().collect { _mseq, index ->
        def target = targets[index] ?: "target_${index + 1}"
        "'${meta.id}_${meta.db_label}_${target}_${index + 1}.mseq'"
    }.join(' ')
    """
    ${target_symlinks}
    mseq_files=( ${input_files} )

    for index in "\${!mseq_files[@]}"
    do
        if [ "\${index}" -eq 0 ]
        then
            cat "\${mseq_files[\${index}]}" > ${merged_mseq}
        else
            # MAPseq writes a header to every file. Keep the first header only.
            awk '!/^#/' "\${mseq_files[\${index}]}" >> ${merged_mseq}
        fi
    done
    """
}
