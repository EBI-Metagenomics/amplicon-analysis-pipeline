
process DADA2 {
    // Run DADA2 pipeline including read-tracking
    tag "$meta.id"
    label "dada2_resources"
    container 'quay.io/microbiome-informatics/dada2:v1'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*map.txt"), path("*asvs.fasta"), path("*_filt.fastq.gz"), optional: true, emit: dada2_out
    tuple val(meta), path("*_dada2_stats.tsv")                                     , optional: true, emit: dada2_stats
    tuple val(meta), path("*_dada2_output.txt")                                    , optional: true, emit: dada2_stdout
    tuple val(meta), path("*_dada2_errors.txt")                                    , optional: true, emit: dada2_errors
    tuple val(meta), path("*_dada2_truncation_points.txt")                         , optional: true, emit: dada2_truncation_points
    tuple val(meta), env(stats_fail)                                               , optional: true, emit: dada2_stats_fail
    path "versions.yml"                                                            , emit: versions
    
    script:
    if ( meta.single_end ){
        """
        output_file="${meta.id}_dada2_output.txt"
        error_file="${meta.id}_dada2_errors.txt"
        dada2.R ${meta.id} $reads 1> \$output_file 2> \$error_file

        grep -E "The (forward|reverse) strand truncation point is: " \$output_file > ${meta.id}_dada2_truncation_points.txt

        stats_fail=false
        if [[ -s \$error_file ]] && grep -q "Caught an error" \$error_file; then
            stats_fail=true
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$( R --version | head -1 | cut -d' ' -f3 )
            dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
        END_VERSIONS
        """
    } else {
        """
        output_file="${meta.id}_dada2_output.txt"
        error_file="${meta.id}_dada2_errors.txt"
        dada2.R ${meta.id} ${reads[0]} ${reads[1]} 1> \$output_file 2> \$error_file

        grep -E "The (forward|reverse) strand truncation point is: " \$output_file > ${meta.id}_dada2_truncation_points.txt

        stats_fail=false
        if [[ -s \$error_file ]] && grep -q "Caught an error" \$error_file; then
            stats_fail=true
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$( R --version | head -1 | cut -d' ' -f3 )
            dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
        END_VERSIONS
        """
    }

}
