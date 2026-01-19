/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT EBI-METAGENOMICS MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { READS_QC } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { READS_QC as READS_QC_MERGE } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { DETECT_RNA } from '../subworkflows/ebi-metagenomics/detect_rna/main'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_SSU } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_LSU } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_PR2 } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_UNITE } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'
include { MAPSEQ_OTU_KRONA as MAPSEQ_OTU_KRONA_ITSONEDB } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MASK_FASTA_SWF } from '../subworkflows/local/mask_fasta_swf.nf'
include { ITS_SANITY_CHECKER } from '../modules/local/its_sanity_checker/main'
include { AMP_REGION_INFERENCE } from '../subworkflows/local/amp_region_inference_swf.nf'
include { PRIMER_IDENTIFICATION } from '../subworkflows/local/primer_identification_swf.nf'
include { AUTOMATIC_PRIMER_PREDICTION } from '../subworkflows/local/automatic_primer_prediction.nf'
include { CONCAT_PRIMER_CUTADAPT } from '../subworkflows/local/concat_primer_cutadapt.nf'
include { DADA2_SWF } from '../subworkflows/local/dada2_swf.nf'
include { MAPSEQ_ASV_KRONA as MAPSEQ_ASV_KRONA_SILVA } from '../subworkflows/local/mapseq_asv_krona_swf.nf'
include { MAPSEQ_ASV_KRONA as MAPSEQ_ASV_KRONA_PR2 } from '../subworkflows/local/mapseq_asv_krona_swf.nf'
include { EXTRACT_ASV_READ_COUNTS } from '../modules/local/extract_asv_read_counts/main'
include { EXTRACT_ASVS_LEFT as EXTRACT_ASVS_LEFT_SILVA } from '../modules/local/extract_asvs_left/main'
include { EXTRACT_ASVS_LEFT as EXTRACT_ASVS_LEFT_PR2 } from '../modules/local/extract_asvs_left/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { DOWNLOAD_FROM_FIRE } from '../modules/ebi-metagenomics/downloadfromfire/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MULTIQC as MULTIQC_RUN } from '../modules/nf-core/multiqc/main.nf'
include { MULTIQC as MULTIQC_STUDY } from '../modules/nf-core/multiqc/main.nf'

// Import dada2 input preparation function (it's very big and deserved to be in its own file) //
include { dada2_input_preparation_function } from '../lib/nf/dada2_input_preparation_function.nf'


// Import samplesheetToList from nf-schema //
include { samplesheetToList } from 'plugin/nf-schema'

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

    // Regular taxonomy resolution method //
    ssu_mapseq_krona_tuple = Channel.value(
        [
            file(params.ssu_db_fasta, checkIfExists: true),
            file(params.ssu_db_tax, checkIfExists: true),
            file(params.ssu_db_otu, checkIfExists: true),
            file(params.ssu_db_mscluster, checkIfExists: true),
            params.ssu_label,
        ]
    )
    lsu_mapseq_krona_tuple = Channel.value(
        [
            file(params.lsu_db_fasta, checkIfExists: true),
            file(params.lsu_db_tax, checkIfExists: true),
            file(params.lsu_db_otu, checkIfExists: true),
            file(params.lsu_db_mscluster, checkIfExists: true),
            params.lsu_label,
        ]
    )
    itsonedb_mapseq_krona_tuple = Channel.value(
        [
            file(params.itsone_db_fasta, checkIfExists: true),
            file(params.itsone_db_tax, checkIfExists: true),
            file(params.itsone_db_otu, checkIfExists: true),
            file(params.itsone_db_mscluster, checkIfExists: true),
            params.itsone_label,
        ]
    )
    unite_mapseq_krona_tuple = Channel.value(
        [
            file(params.unite_db_fasta, checkIfExists: true),
            file(params.unite_db_tax, checkIfExists: true),
            file(params.unite_db_otu, checkIfExists: true),
            file(params.unite_db_mscluster, checkIfExists: true),
            params.unite_label,
        ]
    )
    pr2_mapseq_krona_tuple = Channel.value(
        [
            file(params.pr2_db_fasta, checkIfExists: true),
            file(params.pr2_db_tax, checkIfExists: true),
            file(params.pr2_db_otu, checkIfExists: true),
            file(params.pr2_db_mscluster, checkIfExists: true),
            params.pr2_label,
        ]
    )

    // Regular ASV resolution method //
    dada2_krona_silva_tuple = tuple(
        file(params.ssu_db_fasta, checkIfExists: true),
        file(params.ssu_db_tax, checkIfExists: true),
        file(params.ssu_db_otu, checkIfExists: true),
        file(params.ssu_db_mscluster, checkIfExists: true),
        params.dada2_silva_label,
    )
    dada2_krona_pr2_tuple = tuple(
        file(params.pr2_db_fasta, checkIfExists: true),
        file(params.pr2_db_tax, checkIfExists: true),
        file(params.pr2_db_otu, checkIfExists: true),
        file(params.pr2_db_mscluster, checkIfExists: true),
        params.dada2_pr2_label,
    )


    // Initialise standard primer library for PIMENTO if user-given//
    // If there are no primers provided, it will fallback to use the default PIMENTO standard primer library
    std_primer_library = []

    if (params.std_primer_library) {
        std_primer_library = file(params.std_primer_library, type: 'dir', checkIfExists: true)
    }

    // Read input samplesheet and validate it using schema_input.json //
    samplesheet = Channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))

    ch_versions = Channel.empty()

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

    // Sanity checking and quality control of reads //
    READS_QC_MERGE(
        true,
        ch_input,
        "",
        true,
    )
    ch_versions = ch_versions.mix(READS_QC_MERGE.out.versions)

    // Run it again without merging to keep PE files unmerged for primer trimming+DADA2 //
    READS_QC(
        false,
        ch_input,
        "",
        false,
    )
    ch_versions = ch_versions.mix(READS_QC.out.versions)

    // Removes reads that passed sanity checks but are empty after QC with fastp //
    READS_QC_MERGE.out.reads_fasta
        .branch { _meta, reads ->
            qc_pass: reads.countFasta() > 0
            qc_empty: reads.countFasta() == 0
        }
        .set { extended_reads_qc }

    // rRNA extraction subworkflow to find rRNA reads for SSU+LSU //
    DETECT_RNA(
        extended_reads_qc.qc_pass,
        file(params.rrnas_rfam_covariance_model, checkIfExists: true),
        file(params.rrnas_rfam_claninfo, checkIfExists: true),
        "cmsearch",
        true,
        false,
    )
    ch_versions = ch_versions.mix(DETECT_RNA.out.versions)

    // Masking subworkflow to find rRNA reads for ITS //
    MASK_FASTA_SWF(
        extended_reads_qc.qc_pass,
        DETECT_RNA.out.concat_ssu_lsu_coords,
    )
    ch_versions = ch_versions.mix(MASK_FASTA_SWF.out.versions)

    // Next five subworkflow calls are MAPseq annotation + Krona generation for SSU+LSU+ITS //
    MAPSEQ_OTU_KRONA_SSU(
        DETECT_RNA.out.ssu_fasta,
        ssu_mapseq_krona_tuple,
    )
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA_SSU.out.versions)

    MAPSEQ_OTU_KRONA_PR2(
        DETECT_RNA.out.ssu_fasta,
        pr2_mapseq_krona_tuple,
    )
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA_PR2.out.versions)

    MAPSEQ_OTU_KRONA_LSU(
        DETECT_RNA.out.lsu_fasta,
        lsu_mapseq_krona_tuple,
    )
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA_LSU.out.versions)

    MAPSEQ_OTU_KRONA_ITSONEDB(
        MASK_FASTA_SWF.out.masked_out,
        itsonedb_mapseq_krona_tuple,
    )
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA_ITSONEDB.out.versions)

    MAPSEQ_OTU_KRONA_UNITE(
        MASK_FASTA_SWF.out.masked_out,
        unite_mapseq_krona_tuple,
    )
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA_UNITE.out.versions)

    // Sanity check that masked out ITS reads are not just a different marker gene
    its_sanity_check_input = MASK_FASTA_SWF.out.masked_out
        .join(MAPSEQ_OTU_KRONA_ITSONEDB.out.mseq)
        .join(MAPSEQ_OTU_KRONA_UNITE.out.mseq)

    ITS_SANITY_CHECKER(its_sanity_check_input)


    // Infer amplified variable regions for SSU, extract reads for each amplified region if there are more than one //
    AMP_REGION_INFERENCE(
        DETECT_RNA.out.cmsearch_deoverlap_coords,
        READS_QC_MERGE.out.reads_se_and_merged,
    )
    ch_versions = ch_versions.mix(AMP_REGION_INFERENCE.out.versions)

    // Identify whether primers exist or not in reads, separated by different amplified regions if more than one exists in a run //
    PRIMER_IDENTIFICATION(
        AMP_REGION_INFERENCE.out.extracted_var_out,
        std_primer_library,
    )
    ch_versions = ch_versions.mix(PRIMER_IDENTIFICATION.out.versions)

    // Join primer identification flags with reads belonging to each run+amp_region //
    auto_trimming_input = PRIMER_IDENTIFICATION.out.conductor_out.join(AMP_REGION_INFERENCE.out.extracted_var_out, by: [0])

    /* 
    Run subworkflow for automatic primer prediction
    Outputs empty fasta file if no primers, or fasta file containing predicted primers
    */
    AUTOMATIC_PRIMER_PREDICTION(
        auto_trimming_input
    )
    ch_versions = ch_versions.mix(AUTOMATIC_PRIMER_PREDICTION.out.versions)

    // Concatenate the different combinations of stranded std/auto primers for each run+amp_region //
    concat_input = PRIMER_IDENTIFICATION.out.std_primer_out.join(AUTOMATIC_PRIMER_PREDICTION.out.auto_primer_trimming_out, by: [0])

    // Concatenate all primers for for a run, send them to cutadapt with original QCd reads for primer trimming //
    CONCAT_PRIMER_CUTADAPT(
        concat_input,
        READS_QC.out.reads,
    )
    ch_versions = ch_versions.mix(CONCAT_PRIMER_CUTADAPT.out.versions)


    // Run the large dada2 imput preparation function //
    cutadapt_channel = CONCAT_PRIMER_CUTADAPT.out.cutadapt_out.map { meta, reads ->
        [meta.subMap('id', 'single_end'), meta['var_region'], meta['var_regions_size'], reads]
    }

    dada2_input = dada2_input_preparation_function(concat_input, READS_QC.out.reads, cutadapt_channel)
    // Run DADA2 ASV generation //
    DADA2_SWF(
        dada2_input,
        DETECT_RNA.out.cmsearch_deoverlap_coords,
    )
    ch_versions = ch_versions.mix(DADA2_SWF.out.versions)

    def dada2_stats_fail = DADA2_SWF.out.dada2_stats_fail.map { meta, stats_fail ->
        def key = meta.subMap('id', 'single_end')
        return [key, stats_fail]
    }

    // ASV taxonomic assignments + generate Krona plots for each run+amp_region //
    MAPSEQ_ASV_KRONA_SILVA(
        DADA2_SWF.out.dada2_out,
        AMP_REGION_INFERENCE.out.concat_var_regions,
        AMP_REGION_INFERENCE.out.extracted_var_path,
        dada2_krona_silva_tuple,
    )
    ch_versions = ch_versions.mix(MAPSEQ_ASV_KRONA_SILVA.out.versions)

    MAPSEQ_ASV_KRONA_PR2(
        DADA2_SWF.out.dada2_out,
        AMP_REGION_INFERENCE.out.concat_var_regions,
        AMP_REGION_INFERENCE.out.extracted_var_path,
        dada2_krona_pr2_tuple,
    )
    ch_versions = ch_versions.mix(MAPSEQ_ASV_KRONA_PR2.out.versions)

    /*  
    Multiple steps in ASV calling + annotation can result in lost ASVs
    These final modules make sure the set of ASVs being reported in the different outputs
    are consistent i.e. ASVs in read count files, ASV sequences in FASTA files, etc.
    */
    extract_asv_read_counts_input = MAPSEQ_ASV_KRONA_SILVA.out.asv_read_counts.join(MAPSEQ_ASV_KRONA_PR2.out.asv_read_counts)

    EXTRACT_ASV_READ_COUNTS(extract_asv_read_counts_input)
    ch_versions = ch_versions.mix(EXTRACT_ASV_READ_COUNTS.out.versions)

    extract_asvs_input = EXTRACT_ASV_READ_COUNTS.out.asvs_left
        .filter { meta, asvs_left ->
            meta.var_region != "concat"
        }
        .map { meta, asvs_left ->
            def key = groupKey(meta.subMap('id'), meta.var_regions_size)
            [key, asvs_left]
        }
        .groupTuple(by: 0)
        .join(
            DADA2_SWF.out.dada2_out.map { meta, maps, asv_seqs, filt_reads ->
                [['id': meta.id], asv_seqs]
            }
        )

    extract_asvs_input_silva = extract_asvs_input.join(
        MAPSEQ_ASV_KRONA_SILVA.out.asvtaxtable.map { meta, asvtaxtable ->
            [['id': meta.id], asvtaxtable]
        }
    )
    extract_asvs_input_pr2 = extract_asvs_input.join(
        MAPSEQ_ASV_KRONA_PR2.out.asvtaxtable.map { meta, asvtaxtable ->
            [['id': meta.id], asvtaxtable]
        }
    )

    EXTRACT_ASVS_LEFT_SILVA(extract_asvs_input_silva, dada2_krona_silva_tuple[4])
    ch_versions = ch_versions.mix(EXTRACT_ASVS_LEFT_SILVA.out.versions.first())

    EXTRACT_ASVS_LEFT_PR2(extract_asvs_input_pr2, dada2_krona_pr2_tuple[4])
    ch_versions = ch_versions.mix(EXTRACT_ASVS_LEFT_PR2.out.versions.first())

    /*****************************/
    /* MultiQC reports */
    /****************************/

    // Version collating //
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    multiqc_input = CONCAT_PRIMER_CUTADAPT.out.cutadapt_json
        .map { meta, json ->
            [['id': meta.id, 'single_end': meta.single_end], json]
        }
        .join(READS_QC_MERGE.out.fastp_summary_json, remainder: true)
        .join(
            DADA2_SWF.out.dada2_report.map { meta, tsv ->
                [['id': meta.id, 'single_end': meta.single_end], tsv]
            },
            remainder: true
        )
        .map { meta, cutadapt, fastp, dada2 ->
            def final_inputs = [cutadapt, fastp, dada2]
            if (!cutadapt) {
                final_inputs -= cutadapt
            }
            if (!dada2) {
                final_inputs -= dada2
            }

            [meta, final_inputs]
        }

    // MultiQC for individual runs //
    MULTIQC_RUN(
        multiqc_input,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.first(),
        params.multiqc_config,
        [],
        [],
        [],
        [],
    )

    // MultiQC for study !! assuming we do not have multiple studies in one samplesheet !! //
    multiqc_study = multiqc_input
        .flatten()
        .collect()
        .map { item -> item.findAll { !(it instanceof Map) } }
        .map { dataList ->
            [['id': 'study_multiqc_report'], dataList]
        }

    MULTIQC_STUDY(
        multiqc_study,
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml,
        params.multiqc_config,
        [],
        [],
        [],
        [],
    )

    /*****************************/
    /* End of execution reports */
    /****************************/

    // Extract runs that failed SeqFu check //
    READS_QC.out.seqfu_check
        .splitCsv(sep: "\t", elem: 1)
        .filter { meta, seqfu_res ->
            seqfu_res[0] != "OK"
        }
        .map { meta, __ -> "${meta.id},seqfu_fail" }
        .set { seqfu_fails }

    // Extract runs that failed Suffix Header check //
    READS_QC.out.suffix_header_check
        .filter { meta, sfxhd_res ->
            sfxhd_res.countLines() != 0
        }
        .map { meta, __ -> "${meta.id},sfxhd_fail" }
        .set { sfxhd_fails }

    // Extract runs that failed Library Strategy check //
    READS_QC_MERGE.out.amplicon_check
        .filter { meta, strategy ->
            strategy != "AMPLICON"
        }
        .map { meta, __ -> "${meta.id},libstrat_fail" }
        .set { libstrat_fails }

    // Extract runs that had zero reads after fastp //
    extended_reads_qc.qc_empty
        .map { meta, __ -> "${meta.id},no_reads" }
        .set { no_reads_fails }

    // Save all failed runs to file //
    all_failed_runs = seqfu_fails.concat(sfxhd_fails, libstrat_fails, no_reads_fails)
    all_failed_runs.collectFile(name: "qc_failed_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)

    // Extract passed runs, describe whether those passed runs also ASV results //
    DADA2_SWF.out.dada2_report
        .map { meta, dada2_report -> [["id": meta.id, "single_end": meta.single_end], dada2_report] }
        .concat(extended_reads_qc.qc_pass, dada2_stats_fail)
        .groupTuple()
        .map { meta, results ->
            if (results.size() == 3) {
                return "${meta.id},all_results"
            }
            else {
                if (results[1] == "true") {
                    return "${meta.id},dada2_stats_fail"
                }
                else {
                    return "${meta.id},no_asvs"
                }
            }
            error("Unexpected. meta: ${meta}, results: ${results}")
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

            primer_val.each { run_id, ev, met, gene, region, name, strand, sequence ->
                def new_primer = [
                    "name": name,
                    "region": region,
                    "strand": strand,
                    "sequence": sequence,
                    "identification_strategy": name.contains("_auto") ? "auto" : "std",
                ]
                json_map["primers"] << new_primer
            }

            json_map
        }
        .collect()
        .map { collected_json_maps -> def json_content = new groovy.json.JsonBuilder(collected_json_maps).toPrettyString() }
        .collectFile(name: "primer_validation_summary.json", storeDir: "${params.outdir}", newLine: true, cache: false)
}
