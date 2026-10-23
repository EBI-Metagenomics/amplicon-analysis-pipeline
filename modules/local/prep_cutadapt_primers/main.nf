process PREP_CUTADAPT_PRIMERS {
    tag "$meta.id"
    label 'very_light'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "https://depot.galaxyproject.org/singularity/mgnify-pipelines-toolkit:${params.mpt_version}":
        "biocontainers/mgnify-pipelines-toolkit:${params.mpt_version}" }"

    input:
    tuple val(meta), path(forward_primers), path(reverse_primers)

    output:
    tuple val(meta),
        path("${meta.id}_F.fasta"),
        path("${meta.id}_R.fasta"),
        path("${meta.id}_F_RC.fasta"),
        path("${meta.id}_R_RC.fasta"),
        emit: cutadapt_prepared_primers
    path "versions.yml", emit: versions

    script:
    def n = params.bases_allowed_before_primer
    """
    revcomp_primers.py ${forward_primers} ${meta.id}_F_rc_notprep.fasta
    revcomp_primers.py ${reverse_primers} ${meta.id}_R_rc_notprep.fasta

    prep_cutadapt_primers.py --position start --replace-i --n ${n} ${forward_primers} ${meta.id}_F.fasta
    prep_cutadapt_primers.py --position start --replace-i --n ${n} ${reverse_primers} ${meta.id}_R.fasta
    prep_cutadapt_primers.py --position end --replace-i --n ${n} ${meta.id}_F_rc_notprep.fasta ${meta.id}_F_RC.fasta
    prep_cutadapt_primers.py --position end --replace-i --n ${n} ${meta.id}_R_rc_notprep.fasta ${meta.id}_R_RC.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_F.fasta ${meta.id}_R.fasta ${meta.id}_F_RC.fasta ${meta.id}_R_RC.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgnify-pipelines-toolkit: ${params.mpt_version}
    END_VERSIONS
    """
}
