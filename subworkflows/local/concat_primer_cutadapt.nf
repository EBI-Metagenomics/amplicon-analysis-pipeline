
include { CONCAT_PRIMERS          } from '../../modules/local/concat_primers/main.nf'
include { FINAL_CONCAT_PRIMERS    } from '../../modules/local/final_concat_primers/main.nf'
include { SPLIT_PRIMERS_BY_STRAND } from '../../modules/local/split_primers_by_strand.nf'
include { PRIMER_VALIDATION       } from '../../subworkflows/local/primer_validation_swf.nf'

include { CUTADAPT                } from '../../modules/ebi-metagenomics/cutadapt/main.nf'

workflow CONCAT_PRIMER_CUTADAPT {
    
    take:
        concat_input
        reads
    main:

        ch_versions = Channel.empty()

        CONCAT_PRIMERS(
            concat_input
        )
        ch_versions = ch_versions.mix(CONCAT_PRIMERS.out.versions.first())

        // groupKey solution from https://training.nextflow.io/advanced/grouping/#passing-maps-through-processes
        final_concat_primers_input = CONCAT_PRIMERS.out.concat_primers_out
                                    .map{ meta, concat_primers ->
                                        def key = groupKey(meta.subMap('id', 'single_end'), meta['var_regions_size'])
                                        [ key, meta['var_region'], concat_primers ]
                                    }
                                    .groupTuple(by: 0)
                                    .map{ meta, var_region, final_concat_primers ->
                                        [ meta + ["var_region": var_region], final_concat_primers ]
                                    }
        FINAL_CONCAT_PRIMERS(
            final_concat_primers_input
        )
        ch_versions = ch_versions.mix(FINAL_CONCAT_PRIMERS.out.versions.first())

        primer_validation_input = FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
                              .filter{ meta, primers ->
                                primers.size() > 0
                              }

        runs_without_primers = FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
                              .filter{ meta, primers ->
                                primers.size() == 0
                              }
                              .map{ meta, primers ->
                                [ meta, []]
                              }
    // Verify that any identified primers (both std+auto) actually match to regions of the SSU gene (for Bacteria/Archaea/Eukaryotes)
    // Output of this (a .tsv file) will go to CDCH
        PRIMER_VALIDATION(
            primer_validation_input
        )
        ch_versions = ch_versions.mix(PRIMER_VALIDATION.out.versions)

        post_primer_validation_filter = PRIMER_VALIDATION.out.validated_primers
                                        .mix(runs_without_primers)

        cutadapt_input = post_primer_validation_filter
                        .map{ meta, primers ->
                            [ meta.subMap('id', 'single_end'), meta['var_region'], [meta['var_regions_size']], primers ]
                        }
                        .join(reads, by: [0])
                        .map{ meta, var_region, var_regions_size, primers, final_reads -> 
                            [ meta + ["var_region":var_region, "var_regions_size": var_regions_size], final_reads, primers ]
                        }

        CUTADAPT(
            cutadapt_input
        )
        ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    emit:
        final_concat_primers_out = FINAL_CONCAT_PRIMERS.out.final_concat_primers_out
        cutadapt_out             = CUTADAPT.out.reads
        cutadapt_json            = CUTADAPT.out.json
        primer_validation_out    = PRIMER_VALIDATION.out.primer_validation_out
        versions                 = ch_versions

}