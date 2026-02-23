
include { EXTRACT_RRNA_HITS_FROM_UNMERGED_READS } from '../../modules/local/extract_rrna_hits_from_unmerged_reads'
include { REMOVE_AMBIGUOUS_READS } from '../../modules/local/remove_ambiguous_reads/main.nf'
include { DADA2                  } from '../../modules/local/dada2/main.nf'

workflow DADA2_SWF {
    
    take:
        dada2_input
        cmsearch_deoverlap_out


    main:

        ch_versions = channel.empty()

        seqkit_input = dada2_input
                      .map{ meta, reads ->
                        [ meta.subMap('id', 'single_end'), meta['var_region'], meta['var_regions_size'], reads ]
                       }
                      .join(cmsearch_deoverlap_out, by: [0])
                      .map{meta, var_region, var_regions_size, reads, cmsearch_deoverlap_out_ ->
                        [meta + ["var_region": var_region, "var_regions_size": var_regions_size], reads, cmsearch_deoverlap_out_]
                      }

        EXTRACT_RRNA_HITS_FROM_UNMERGED_READS(seqkit_input)

        REMOVE_AMBIGUOUS_READS(
            EXTRACT_RRNA_HITS_FROM_UNMERGED_READS.out.extracted_reads
        )
        ch_versions = ch_versions.mix(REMOVE_AMBIGUOUS_READS.out.versions.first())

        DADA2(
            REMOVE_AMBIGUOUS_READS.out.noambig_out
        )
        ch_versions = ch_versions.mix(DADA2.out.versions.first())

    emit:
        dada2_out        = DADA2.out.dada2_out
        dada2_report     = DADA2.out.dada2_stats
        dada2_errors     = DADA2.out.dada2_errors
        dada2_stats_fail = DADA2.out.dada2_stats_fail
        versions         = ch_versions
    
}