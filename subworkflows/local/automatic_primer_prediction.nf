
include { PIMENTO_GENERATEBCV   } from '../../modules/ebi-metagenomics/pimento/generatebcv/main'
include { FIND_PRIMER_CUTOFFS   } from '../../modules/local/find_primer_cutoffs/main.nf'
include { CHOOSE_PRIMER_CUTOFF } from '../../modules/local/choose_primer_cutoff/main.nf'

workflow AUTOMATIC_PRIMER_PREDICTION {

    // Subworkflow that performs automatic primer trimming on sequencing reads
    // Takes flags for the forward and reverse strands to make sure primers are  
    // predicted only for the required strands
    // In theory this should work on single-end reads as 'forward' and 'reverse'
    // more accurately means 5' or 3' ends of a merged paired-end file
    // Outputs a fasta file that is either empty or that contains the predicted primer sequences (5' to 3')
    
    take:
        auto_trimming_input
    main:

        ch_versions = channel.empty()
        // Use Most Common Prefix (MCP) method to generate curves of base conservation
        PIMENTO_GENERATEBCV(
            auto_trimming_input,
        )
        ch_versions = ch_versions.mix(PIMENTO_GENERATEBCV.out.versions.first())

        // Find inflection points in conservation curves
        FIND_PRIMER_CUTOFFS(
            PIMENTO_GENERATEBCV.out.tsv
        )
        ch_versions = ch_versions.mix(FIND_PRIMER_CUTOFFS.out.versions.first())

        extracted_reads = auto_trimming_input.map{ meta, _fwd_flag, _rev_flag, extracted_reads ->
                                                    [ meta, extracted_reads ]
                                                 }
        // Join fastq channel and the inf_points channel
        assess_inf_input = FIND_PRIMER_CUTOFFS.out.cutoffs
                           .join(extracted_reads, by: [0])

        // Select inflection points most likely to be primer cutoff points
        CHOOSE_PRIMER_CUTOFF(
            assess_inf_input
        )
        ch_versions = ch_versions.mix(CHOOSE_PRIMER_CUTOFF.out.versions.first())

        CHOOSE_PRIMER_CUTOFF.out.auto_primer_out

   emit:
        auto_primer_trimming_out = CHOOSE_PRIMER_CUTOFF.out.auto_primer_out
        versions = ch_versions

}