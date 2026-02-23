
include { STD_PRIMER_FLAG     } from '../../modules/local/std_primer_flag/main.nf'
include { GENERAL_PRIMER_FLAG } from '../../modules/local/general_primer_flag/main.nf'
include { TRIMMING_CONDUCTOR  } from '../../modules/local/trimming_conductor/main.nf'
include { PARSE_CONDUCTOR     } from '../../modules/local/parse_conductor/main.nf'

workflow PRIMER_IDENTIFICATION {
    
    // Subworkflow that attempts to identify the presence of primers in sequencing files (fastq)
    // Checks for the presence of standard primers from a known library in ""./data/standard_primers"
    // Also checks for the presence of primers in general in case a presence is present but is not part of the current library
    // Outputs strand-based flags outlining whether primers were found, and if so, how they should be removed 
    // This can be either removal through matches to known standard primers, or the need to predict the used primer if it's not part of the current library

    take:
        reads_merged
        std_primer_library
    main:

        ch_versions = channel.empty()

        // Standard Library primers
        STD_PRIMER_FLAG(
            reads_merged,
            std_primer_library
        )
        ch_versions = ch_versions.mix(STD_PRIMER_FLAG.out.versions.first())

        // Primers in general
        GENERAL_PRIMER_FLAG(
            reads_merged
        )
        ch_versions = ch_versions.mix(GENERAL_PRIMER_FLAG.out.versions.first())

        // Combining std and general primer outputs and parsing them to guide 
        // samples through automatic primer identification and trimming by cutadapt
        comb_flags = GENERAL_PRIMER_FLAG.out.general_primer_out
                     .join(STD_PRIMER_FLAG.out.std_primer_out, by: [0])
        
        // Generate a flags file with vals of 'none/std/auto' for both fwd and rev
        TRIMMING_CONDUCTOR(
            comb_flags
        )
        ch_versions = ch_versions.mix(TRIMMING_CONDUCTOR.out.versions.first())
        
        // Parse flags file into env variables
        PARSE_CONDUCTOR(
            TRIMMING_CONDUCTOR.out.trimming_conductor_out
        )
        ch_versions = ch_versions.mix(PARSE_CONDUCTOR.out.versions.first())

    emit:
        conductor_out = PARSE_CONDUCTOR.out.conductor_out
        std_primer_out = STD_PRIMER_FLAG.out.std_primer_out
        versions = ch_versions
    
}