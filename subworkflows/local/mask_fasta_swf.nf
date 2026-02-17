
include { FORMAT_BEDFILE     } from '../../modules/local/format_bedfile/main'
include { BEDTOOLS_MASKFASTA } from '../../modules/nf-core/bedtools/maskfasta/main'
include { FILTER_MASKED_N } from '../../modules/local/filter_masked_n/main'

workflow MASK_FASTA_SWF {
    
    take:
        reads_fasta
        concat_ssu_lsu_coords
    main:

        ch_versions = Channel.empty()

        FORMAT_BEDFILE(
            concat_ssu_lsu_coords,
        )
        ch_versions = ch_versions.mix(FORMAT_BEDFILE.out.versions.first())

        bedtools_input = reads_fasta
                         .join(FORMAT_BEDFILE.out.format_bedfile_out)

        BEDTOOLS_MASKFASTA(
            bedtools_input
        )
        ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions.first())

        FILTER_MASKED_N(BEDTOOLS_MASKFASTA.out.fasta)
        ch_versions = ch_versions.mix(FILTER_MASKED_N.out.versions.first())


    emit:
        masked_out = FILTER_MASKED_N.out.filtered_its_fasta
        versions = ch_versions
}