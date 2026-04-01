
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
    tuple val(meta), path("*_dada2_errors.txt")                                    , optional: true, emit: dada2_errors
    tuple val(meta), path("*_dada2_truncation_points.tsv")                         , optional: true, emit: dada2_truncation_points
    tuple val(meta), env(stats_fail)                                               , optional: true, emit: dada2_stats_fail
    path "versions.yml"                                                            , emit: versions
    
    script:
    if ( meta.single_end ){
        """
        stdout_file="${meta.id}_dada2_stdout.txt"
        stderr_file="${meta.id}_dada2_stderr.txt"
        error_file="${meta.id}_dada2_errors.txt"
        dada2.R ${meta.id} $reads 1> \$stdout_file 2> \$stderr_file

        fwd_trunc=\$(grep -E "The forward strand truncation point is: " \$stdout_file | sed -E 's/^\\[[0-9]+\\][^0-9]*([0-9]+)[^0-9]*\$/\\1/')
        rev_trunc=\$(grep -E "The reverse strand truncation point is: " \$stdout_file | sed -E 's/^\\[[0-9]+\\][^0-9]*([0-9]+)[^0-9]*\$/\\1/')
        
        trunc_fp="${meta.id}_dada2_truncation_points.tsv"
        echo "pair\tposition" > \$trunc_fp
        if [ -n \$fwd_trunc ]; then
            echo "forward\t\$fwd_trunc" >> \$trunc_fp
        fi
        if [ -n \$rev_trunc ]; then
            echo "reverse\t\$rev_trunc" >> \$trunc_fp
        fi

        stats_fail=false
        if [[ -s \$stderr_file ]] && grep -q "Caught an error" \$stderr_file; then
            stats_fail=true
            mv \$stderr_file \$error_file
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$( R --version | head -1 | cut -d' ' -f3 )
            dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
        END_VERSIONS
        """
    } else {
        """
        stdout_file="${meta.id}_dada2_stdout.txt"
        stderr_file="${meta.id}_dada2_stderr.txt"
        error_file="${meta.id}_dada2_errors.txt"
        dada2.R ${meta.id} ${reads[0]} ${reads[1]} 1> \$stdout_file 2> \$stderr_file

        fwd_trunc=\$(grep -E "The forward strand truncation point is: " \$stdout_file | sed -E 's/^\\[[0-9]+\\][^0-9]*([0-9]+)[^0-9]*\$/\\1/')
        rev_trunc=\$(grep -E "The reverse strand truncation point is: " \$stdout_file | sed -E 's/^\\[[0-9]+\\][^0-9]*([0-9]+)[^0-9]*\$/\\1/')
        
        trunc_fp="${meta.id}_dada2_truncation_points.tsv"
        echo "pair\tposition" > \$trunc_fp
        if [ -n \$fwd_trunc ]; then
            echo "forward\t\$fwd_trunc" >> \$trunc_fp
        fi
        if [ -n \$rev_trunc ]; then
            echo "reverse\t\$rev_trunc" >> \$trunc_fp
        fi

        stats_fail=false
        if [[ -s \$stderr_file ]] && grep -q "Caught an error" \$stderr_file; then
            stats_fail=true
            mv \$stderr_file \$error_file
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: \$( R --version | head -1 | cut -d' ' -f3 )
            dada2: \$( R -e "suppressMessages(library(dada2));packageDescription('dada2')" | grep 'Version' | cut -d' ' -f2 )
        END_VERSIONS
        """
    }

}
