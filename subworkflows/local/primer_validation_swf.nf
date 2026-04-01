
include { PERMUTE_PRIMERS                                        } from '../../modules/local/permute_primers/main.nf'
include { INFERNAL_CMSEARCH as PRIMER_VALIDATION_SEARCH          } from '../../modules/ebi-metagenomics/infernal/cmsearch/main.nf'
include { CMSEARCHTBLOUTDEOVERLAP as PRIMER_VALIDATION_DEOVERLAP } from '../../modules/ebi-metagenomics/cmsearchtbloutdeoverlap/main.nf'
include { PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS                 } from '../../modules/local/primer_validation_classify_var_regions/main.nf'

workflow PRIMER_VALIDATION {
    
    take:
        primer_validation_input

    main:

        ch_versions = channel.empty()

        PERMUTE_PRIMERS(primer_validation_input)
        ch_versions = ch_versions.mix(PERMUTE_PRIMERS.out.versions.first())

        PRIMER_VALIDATION_SEARCH(
            PERMUTE_PRIMERS.out.permuted_primers,
            file(params.rrnas_rfam_covariance_model, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(PRIMER_VALIDATION_SEARCH.out.versions.first())
        
        PRIMER_VALIDATION_DEOVERLAP(
            PRIMER_VALIDATION_SEARCH.out.cmsearch_tbl,
            file(params.rrnas_rfam_claninfo, checkIfExists: true)
        )
        ch_versions = ch_versions.mix(PRIMER_VALIDATION_DEOVERLAP.out.versions.first())

        primer_validation_classify_var_regions_input = PRIMER_VALIDATION_DEOVERLAP.out.cmsearch_tblout_deoverlapped
                                                       .map{ it -> tuple(["id":it[0].id, "single_end":it[0].single_end], it[0].var_region, it[1]) }
                                                       .join(PERMUTE_PRIMERS.out.permuted_primers.map { meta, permuted_primers -> {
                                                            tuple(["id": meta.id, "single_end": meta.single_end], meta.var_region, permuted_primers)
                                                            }
                                                        }, by: 0)
                                                       .map{ it -> tuple( it[0] + ["var_region":it[1]], it[2], it[4]) }
        
        PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS(
            primer_validation_classify_var_regions_input
        )
        ch_versions = ch_versions.mix(PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS.out.versions.first())

    emit:
        primer_validation_out = PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS.out.primer_validation_out
        validated_primers     = PRIMER_VALIDATION_CLASSIFY_VAR_REGIONS.out.validated_primers
        versions              = ch_versions
    
}