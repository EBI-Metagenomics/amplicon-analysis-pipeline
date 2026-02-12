/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT EBI-METAGENOMICS MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READS_QC                         } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { READS_QC as READS_QC_MERGE       } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { DETECT_RNA                       } from '../subworkflows/ebi-metagenomics/detect_rna/main'
include { MAPSEQ_OTU_KRONA                 } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BBMAP_REFORMAT_STANDARDISE       } from '../modules/local/bbmap/reformat_standardise/main'
include { MASK_FASTA_SWF                   } from '../subworkflows/local/mask_fasta_swf.nf'
include { AMP_REGION_INFERENCE             } from '../subworkflows/local/amp_region_inference_swf.nf'
include { PRIMER_IDENTIFICATION            } from '../subworkflows/local/primer_identification_swf.nf'
include { AUTOMATIC_PRIMER_PREDICTION      } from '../subworkflows/local/automatic_primer_prediction.nf'
include { CONCAT_PRIMER_CUTADAPT           } from '../subworkflows/local/concat_primer_cutadapt.nf'
include { DADA2_SWF                        } from '../subworkflows/local/dada2_swf.nf'
include { MAPSEQ_ASV_KRONA                 } from '../subworkflows/local/mapseq_asv_krona_swf.nf'
include { EXTRACT_ASV_READ_COUNTS          } from '../modules/local/extract_asv_read_counts/main'
include { EXTRACT_ASVS_LEFT                } from '../modules/local/extract_asvs_left/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_FROM_FIRE               } from '../modules/ebi-metagenomics/downloadfromfire/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS      } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC as MULTIQC_RUN           } from '../modules/nf-core/multiqc/main.nf'
include { MULTIQC as MULTIQC_STUDY         } from '../modules/nf-core/multiqc/main.nf'

// Import dada2 input preparation function (it's very big and deserved to be in its own file) //
include { dada2_input_preparation_function } from '../lib/nf/dada2_input_preparation_function.nf'


// Import samplesheetToList from nf-schema //
include { samplesheetToList                } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow AMPLICON_PIPELINE {

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INITIALISE REFERENCE DATABASE INPUT TUPLES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // Parse config to get amplicon reference databases
    mapseq_dbs_in = channel
        .from(
            params.mapseq_databases.collect { k, v ->
                if (v instanceof Map) {
                    if (v.containsKey('label')) {
                        return [[id: k], v]
                    }
                }
            }
        )
        .filter { it -> it }
        .map { meta, fields -> 
            def asv_label = fields.run_asv ? ['asv_label': fields.asv_label] : []
            def extra_meta = ['label': fields.label, 'asv': fields.run_asv, 'otu': fields.run_otu]
            [
                meta + extra_meta + asv_label, 
                tuple(
                    file(fields.fasta), 
                    file(fields.tax), 
                    file(fields.otu), 
                    file(fields.mscluster), 
                    fields.label
                )
            ]
        }
    mapseq_dbs_in.view{ it -> "mapseq_dbs_in - ${it}"}

    // Initialise standard primer library for PIMENTO if user-given//
    // If there are no primers provided, it will fallback to use the default PIMENTO standard primer library
    std_primer_library = []

    if (params.std_primer_library){
        std_primer_library = file(params.std_primer_library, type: 'dir', checkIfExists: true)
    }

    // Read input samplesheet and validate it using schema_input.json //
    samplesheet = channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))

    ch_versions = channel.empty()

    // Organise input tuple channel //
    groupReads = { meta, fq1, fq2 ->
        if (fq2 == []) {
            return tuple(meta, [fq1])
        }
        else {
            return tuple(meta, [fq1, fq2])
        }
    }

    ch_input = samplesheet.map(groupReads)

    if (params.use_fire_download) {
        /*
         * Internally we need to bypass Nextflow S3 integration until https://github.com/nextflow-io/nextflow/issues/4873 is fixed
         * The EBI parameter is needed as this only works on EBI network, FIRE is not accessible otherwise
        */
        DOWNLOAD_FROM_FIRE(
            ch_input
        )

        ch_versions = ch_versions.mix(DOWNLOAD_FROM_FIRE.out.versions.first())
        ch_input = DOWNLOAD_FROM_FIRE.out.downloaded_files
    }

    // De-interleave interleaved paired-end reads
    BBMAP_REFORMAT_STANDARDISE(ch_input, 'fastq.gz')
    ch_input = BBMAP_REFORMAT_STANDARDISE.out.reformated

    // Sanity checking and quality control of reads //
    READS_QC_MERGE(
        true, // check if amplicon
        ch_input,
        "",  // don't discard trimmed reads
        true // merge
    )
    ch_versions = ch_versions.mix(READS_QC_MERGE.out.versions)

    // Run it again without merging to keep PE files unmerged for primer trimming+DADA2 //
    READS_QC(
        false, // check if amplicon
        ch_input,
        "",  // don't discard trimmed reads
        false // merge
    )
    ch_versions = ch_versions.mix(READS_QC.out.versions)

    // Removes reads that passed sanity checks but are empty after QC with fastp //
    READS_QC_MERGE.out.reads_fasta
        .branch{ _meta, reads ->
            qc_pass: reads.countFasta() > 0
            qc_empty: reads.countFasta() == 0
        }
        .set { extended_reads_qc }

    // rRNA extraction subworkflow to find rRNA reads for SSU+LSU //
    DETECT_RNA(
        extended_reads_qc.qc_pass,
        file( params.rrnas_rfam_covariance_model, checkIfExists: true ),
        file( params.rrnas_rfam_claninfo, checkIfExists: true ),
        "cmsearch",
        true,
        false
    )
    ch_versions = ch_versions.mix(DETECT_RNA.out.versions)

    // Masking subworkflow to find rRNA reads for ITS //
    MASK_FASTA_SWF(
        extended_reads_qc.qc_pass,
        DETECT_RNA.out.concat_ssu_lsu_coords
    )
    ch_versions = ch_versions.mix(MASK_FASTA_SWF.out.versions)

    
    // CHANGE HERE
    // Next five subworkflow calls are MAPseq annotation + Krona generation for SSU+LSU+ITS //

    mapseq_otu_dbs_in = mapseq_dbs_in.filter{ meta, _db -> meta.otu }
    MAPSEQ_OTU_KRONA(DETECT_RNA.out.ssu_fasta, mapseq_otu_dbs_in)
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA.out.versions)


    if (!params.skip_asv) {
        // Infer amplified variable regions for SSU, extract reads for each amplified region if there are more than one //
        AMP_REGION_INFERENCE(
            DETECT_RNA.out.cmsearch_deoverlap_coords,
            READS_QC_MERGE.out.reads_se_and_merged
        )
        ch_versions = ch_versions.mix(AMP_REGION_INFERENCE.out.versions)

        // Identify whether primers exist or not in reads, separated by different amplified regions if more than one exists in a run //
        PRIMER_IDENTIFICATION(
            AMP_REGION_INFERENCE.out.extracted_var_out,
            std_primer_library
        )
        ch_versions = ch_versions.mix(PRIMER_IDENTIFICATION.out.versions)

        // Join primer identification flags with reads belonging to each run+amp_region //
        auto_trimming_input = PRIMER_IDENTIFICATION.out.conductor_out
                              .join(AMP_REGION_INFERENCE.out.extracted_var_out, by: [0])

        /* 
        Run subworkflow for automatic primer prediction
        Outputs empty fasta file if no primers, or fasta file containing predicted primers
        */
        AUTOMATIC_PRIMER_PREDICTION(
            auto_trimming_input
        )
        ch_versions = ch_versions.mix(AUTOMATIC_PRIMER_PREDICTION.out.versions)

        // Concatenate the different combinations of stranded std/auto primers for each run+amp_region //
        concat_input = PRIMER_IDENTIFICATION.out.std_primer_out
                       .join(AUTOMATIC_PRIMER_PREDICTION.out.auto_primer_trimming_out, by: [0])

        // Concatenate all primers for for a run, send them to cutadapt with original QCd reads for primer trimming //
        CONCAT_PRIMER_CUTADAPT(
            concat_input,
            READS_QC.out.reads
        )
        ch_versions = ch_versions.mix(CONCAT_PRIMER_CUTADAPT.out.versions)


        // Run the large dada2 input preparation function //
        cutadapt_channel = CONCAT_PRIMER_CUTADAPT.out.cutadapt_out
            .map { meta, reads -> 
                [ meta.subMap('id', 'single_end'), meta['var_region'], meta['var_regions_size'], reads ]
            }

        dada2_input = dada2_input_preparation_function(concat_input, READS_QC.out.reads, cutadapt_channel)
        // Run DADA2 ASV generation //
        DADA2_SWF(
            dada2_input,
            DETECT_RNA.out.cmsearch_deoverlap_coords
        )
        ch_versions = ch_versions.mix(DADA2_SWF.out.versions)


        // CHANGE FROM HERE
        // ASV taxonomic assignments + generate Krona plots for each run+amp_region //
        mapseq_asv_dbs_in = mapseq_dbs_in.filter{ meta, _db -> meta.asv }
        MAPSEQ_ASV_KRONA(
            DADA2_SWF.out.dada2_out,
            AMP_REGION_INFERENCE.out.concat_var_regions,
            AMP_REGION_INFERENCE.out.extracted_var_path,
            mapseq_asv_dbs_in,
        )
        ch_versions = ch_versions.mix(MAPSEQ_ASV_KRONA.out.versions)

        /*  
        Multiple steps in ASV calling + annotation can result in lost ASVs
        These final modules make sure the set of ASVs being reported in the different outputs
        are consistent i.e. ASVs in read count files, ASV sequences in FASTA files, etc.
        */
        extract_asv_read_counts_input = MAPSEQ_ASV_KRONA.out.asv_read_counts
            .map{ meta, counts -> [meta.subMap('id', 'var_region', 'var_regions_size', 'asv_label'), counts] }
            .groupTuple()
        EXTRACT_ASV_READ_COUNTS(extract_asv_read_counts_input)
        ch_versions = ch_versions.mix(EXTRACT_ASV_READ_COUNTS.out.versions)
        
        extract_asvs_input = EXTRACT_ASV_READ_COUNTS.out.asvs_left
            .filter { meta, _asvs_left -> meta.var_region != "concat" }
            .map{ meta, asvs_left ->
                def renamed_meta = ['id': meta.id, 'asv_label': meta.asv_label]
                def key = groupKey(renamed_meta, meta.var_regions_size)
                return [ key, asvs_left ]
            }
            .groupTuple()
            .map{ meta, asvs_left -> [meta.subMap('id'), meta, asvs_left] }
            .combine(
                DADA2_SWF.out.dada2_out
                    .map { meta, _maps, asv_seqs, _filt_reads ->
                           [meta.subMap('id'), meta, asv_seqs] },
                by: 0
            )
            .map{ _meta_id, meta, asvs_left, _dada2_meta, asv_seqs -> 
                  [meta, asvs_left, asv_seqs] }
            .join(MAPSEQ_ASV_KRONA.out.asvtaxtable
                .map{ meta, asvtaxtable ->
                      [meta.subMap('id', 'asv_label'), asvtaxtable] }
            )
            .map{ meta, asvs_left, asv_seqs, asvtaxtable -> 
                  [meta, asvs_left, asv_seqs, asvtaxtable, meta.asv_label] }
        EXTRACT_ASVS_LEFT(extract_asvs_input)
        ch_versions = ch_versions.mix(EXTRACT_ASVS_LEFT.out.versions.first())


        /*****************************/
        /* End of execution reports */
        /****************************/

        def dada2_stats_fail = DADA2_SWF.out.dada2_stats_fail
            .map { meta, stats_fail ->
                def key = meta.subMap('id', 'single_end')
                return [key, stats_fail]
            }

        // Extract passed runs, describe whether those passed runs also ASV results //
        DADA2_SWF.out.dada2_report
            .map { meta, dada2_report -> [ ["id": meta.id, "single_end": meta.single_end], dada2_report ] }
            .concat(extended_reads_qc.qc_pass, dada2_stats_fail)
            .groupTuple()
            .map { meta, results ->
                if ( results.size() == 3 ) {
                    return "${meta.id},all_results"
                }
                else {
                    if (results[1] == "true"){
                        return "${meta.id},dada2_stats_fail"
                    } else {
                        return "${meta.id},no_asvs"
                    }
                }
                error "Unexpected. meta: ${meta}, results: ${results}"
            }
            .set { final_passed_runs }

        // Save all passed runs to file //
        final_passed_runs
            .collectFile(name: "qc_passed_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)
            .set { passed_runs_path }

        // Summarise primer validation information into study-wide JSON file //
        CONCAT_PRIMER_CUTADAPT.out.primer_validation_out
            .splitCsv(sep: "\t", elem: 1, skip: 1)
            .groupTuple()
            .map { meta, primer_val ->

                def json_map = ["id": "${meta.id}", "primers": []]

                primer_val.each { _run_id, _ev, _met, _gene, region, name, strand, sequence ->
                    def new_primer = [
                        "name": name,
                        "region": region,
                        "strand": strand,
                        "sequence": sequence,
                        "identification_strategy": name.contains("_auto") ? "auto" : "std"
                    ]
                    json_map["primers"] << new_primer
                }

                return json_map
             }
            .collect()
            .map { collected_json_maps -> 
                def json_content = new groovy.json.JsonBuilder(collected_json_maps).toPrettyString() }
            .collectFile(
                name: "primer_validation_summary.json", 
                storeDir: "${params.outdir}", 
                newLine: true, 
                cache: false
            )
    } 

    /*****************************/
    /* MultiQC reports */
    /****************************/

    // Version collating //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    if (params.skip_asv) {
        multiqc_input = READS_QC_MERGE.out.fastp_summary_json
                        .map{ meta, fastp ->
                                def final_inputs = [fastp]
                                [meta, final_inputs]
                            }
    } else {
        multiqc_input = CONCAT_PRIMER_CUTADAPT.out.cutadapt_json.map{ meta, json ->
                            [['id':meta.id, 'single_end':meta.single_end], json]
                        }
                        .join(READS_QC_MERGE.out.fastp_summary_json, remainder:true)
                        .join(DADA2_SWF.out.dada2_report.map{ meta, tsv ->
                            [['id':meta.id, 'single_end':meta.single_end], tsv]}, remainder:true)
                        .map{ meta, cutadapt, fastp, dada2 ->
                                def final_inputs = [cutadapt, fastp, dada2]
                                if (!cutadapt){
                                    final_inputs -= cutadapt
                                }
                                if (!dada2){
                                    final_inputs -= dada2
                                }

                                [meta, final_inputs]
                            }
    }

    // MultiQC for individual runs //
    MULTIQC_RUN(multiqc_input,
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.first(),
            params.multiqc_config,
            [],
            [],
            [],
            []
            )

    // MultiQC for study !! assuming we do not have multiple studies in one samplesheet !! //
    multiqc_study = multiqc_input.flatten().collect()
        .map{ item ->item.findAll { it -> !(it instanceof Map) }}
        .map { dataList ->
            [['id': 'study_multiqc_report'], dataList ]
        }

    MULTIQC_STUDY(multiqc_study,
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml,
            params.multiqc_config,
            [],
            [],
            [],
            []
            )

    /*****************************/
    /* End of execution reports */
    /****************************/

    // Extract runs that failed SeqFu check //
    READS_QC.out.seqfu_check
        .splitCsv(sep: "\t", elem: 1)
        .filter { _meta, seqfu_res ->
            seqfu_res[0] != "OK"
        }
        .map { meta, __ -> "${meta.id},seqfu_fail" }
        .set { seqfu_fails }

    // Extract runs that failed Suffix Header check //
    READS_QC.out.suffix_header_check
        .filter { _meta, sfxhd_res ->
            sfxhd_res.countLines() != 0
        }
        .map { meta, __ -> "${meta.id},sfxhd_fail"  }
        .set { sfxhd_fails }

    // Extract runs that failed Library Strategy check //
    READS_QC_MERGE.out.amplicon_check
        .filter { _meta, strategy ->
            strategy != "AMPLICON"
        }
        .map { meta, __ -> "${meta.id},libstrat_fail" }
        .set { libstrat_fails }

    // Extract runs that had zero reads after fastp //
    extended_reads_qc.qc_empty.map { meta, __ -> "${meta.id},no_reads"  }
        .set { no_reads_fails }

    // Save all failed runs to file //
    all_failed_runs = seqfu_fails.concat( sfxhd_fails, libstrat_fails, no_reads_fails )
    all_failed_runs.collectFile(
        name: "qc_failed_runs.csv", 
        storeDir: "${params.outdir}", 
        newLine: true, 
        cache: false
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
