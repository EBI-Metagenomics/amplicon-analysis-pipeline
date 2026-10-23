/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT EBI-METAGENOMICS MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BBMAP_REFORMAT_STANDARDISE       } from '../modules/ebi-metagenomics/bbmap/reformat_standardise/main'
include { READS_QC                         } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { READS_QC as READS_QC_MERGE       } from '../subworkflows/ebi-metagenomics/reads_qc/main.nf'
include { DETECT_RNA                       } from '../subworkflows/ebi-metagenomics/detect_rna/main'
include { MAPSEQ_OTU_KRONA                 } from '../subworkflows/ebi-metagenomics/mapseq_otu_krona/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MASK_FASTA_SWF                   } from '../subworkflows/local/mask_fasta_swf.nf'
include { AMP_REGION_INFERENCE             } from '../subworkflows/local/amp_region_inference_swf.nf'
include { PRIMER_IDENTIFICATION            } from '../subworkflows/local/primer_identification_swf.nf'
include { AUTOMATIC_PRIMER_PREDICTION      } from '../subworkflows/local/automatic_primer_prediction.nf'
include { CONCAT_PRIMER_CUTADAPT           } from '../subworkflows/local/concat_primer_cutadapt.nf'
include { DADA2_SWF                        } from '../subworkflows/local/dada2_swf.nf'
include { MAPSEQ_ASV_KRONA                 } from '../subworkflows/local/mapseq_asv_krona_swf.nf'
include { MAKE_ASV_COUNT_TABLES            } from '../modules/local/make_asv_count_tables/main'
include { ITS_SANITY_CHECKER               } from '../modules/local/its_sanity_checker/main'
include { FILTER_MAPSEQ_BY_MASK            } from '../modules/local/filter_mapseq_by_mask/main'
include { MERGE_MAPSEQ                     } from '../modules/local/merge_mapseq/main'
include { PUBLISH_OTU_RESULTS              } from '../modules/local/publish_otu_results/main'
include { PUBLISH_MASKED_ITS               } from '../modules/local/publish_masked_its/main'

include { MAPSEQ2BIOM as MAPSEQ2BIOM_MERGED } from '../modules/ebi-metagenomics/mapseq2biom/main'
include { KRONA_KTIMPORTTEXT as KRONA_KTIMPORTTEXT_MERGED } from '../modules/ebi-metagenomics/krona/ktimporttext/main'

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
            if (fields.run_asv) {
                if (!fields.containsKey('asv_label')) {
                    throw new IllegalArgumentException("ASV database '${meta.id}' has run_asv=true but no asv_label")
                }
                if (!fields.containsKey('tax_ranks')) {
                    throw new IllegalArgumentException("ASV database '${meta.id}' has run_asv=true but no tax_ranks")
                }
                if (!(fields.tax_ranks in ['DADA2-SILVA', 'DADA2-PR2'])) {
                    throw new IllegalArgumentException("ASV database '${meta.id}' has unsupported tax_ranks '${fields.tax_ranks}'. Expected DADA2-SILVA or DADA2-PR2")
                }
            }
            def asv_meta = fields.run_asv ? ['asv_label': fields.asv_label, 'tax_ranks': fields.tax_ranks] : []
            def extra_meta = [
                'target': fields.target,
                'db_label': fields.label,
                'run_asv': fields.run_asv,
                'run_otu': fields.run_otu,
            ]

            [
                meta + extra_meta + asv_meta, 
                tuple(
                    file(fields.fasta), 
                    file(fields.tax), 
                    file(fields.otu), 
                    file(fields.mscluster), 
                    fields.label
                )
            ]
        }

    // Initialise standard primer library for PIMENTO if user-given//
    // If there are no primers provided, it will fallback to use the default PIMENTO standard primer library
    std_primer_library = []

    if (params.std_primer_library) {
        std_primer_library = file(params.std_primer_library, type: 'dir', checkIfExists: true)
    }

    // Read input samplesheet and validate it using schema_input.json //
    samplesheet = channel.fromList(samplesheetToList(params.input, "./assets/schema_input.json"))

    ch_versions = channel.empty()

    // Organise input tuple channel //
    def groupReads = { meta, fq1, fq2 ->
        def single_file = (fq2 == [])
        meta['interleaved'] = (!meta.single_end) && single_file
        if (single_file) {
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
    
    // Standardise headers and de-interleave as needed
    if (!params.skip_standardise) {
        standardise_input = ch_input.multiMap{
            meta, reads ->
            reads: [meta, reads]
            interleaved: meta.interleaved
        }
        BBMAP_REFORMAT_STANDARDISE(
            standardise_input.reads, 
            standardise_input.interleaved, 
            'fastq.gz'
        )
        ch_versions = ch_versions.mix(BBMAP_REFORMAT_STANDARDISE.out.versions)
        ch_input = BBMAP_REFORMAT_STANDARDISE.out.reformated
    }
    
    // Sanity checking and quality control of reads //
    READS_QC_MERGE(
        true,       // check if amplicon
        ch_input,
        "",         // don't discard trimmed reads
        true,       // merge
    )
    ch_versions = ch_versions.mix(READS_QC_MERGE.out.versions)

    // Run it again without merging to keep PE files unmerged for primer trimming+DADA2 //
    READS_QC(
        false,      // don't check if amplicon
        ch_input,
        "",         // don't discard trimmed reads
        false,      // don't merge
    )
    ch_versions = ch_versions.mix(READS_QC.out.versions)

    // Removes reads that passed sanity checks but are empty after QC with fastp //
    extended_reads_qc = READS_QC_MERGE.out.reads_fasta
        .branch { _meta, reads ->
            qc_pass: reads.countFasta() > 0
            qc_empty: reads.countFasta() == 0
        }

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
        DETECT_RNA.out.all_identified_coords,
    )
    ch_versions = ch_versions.mix(MASK_FASTA_SWF.out.versions)

    
    // The next subworkflow is MAPseq annotation + Krona generation for SSU+LSU+ITS //

    // We launch it on all generated fasta files to which we associate one or two target databases:
    // ssu_fasta for SSU and potential ITS (as 16S and 18S are followed by ITS)
    // lsu_fasta for LSU and potential ITS (as 23S and 25S/28S are preceded by ITS)
    // five_eightS_fasta for potential ITS (as it's inbetween ITS1 and ITS2)

    ch_ssu_fasta = DETECT_RNA.out.ssu_fasta
        .map { meta, fasta ->
            tuple(meta + [target: 'SSU'], fasta)
        }

    ch_lsu_fasta = DETECT_RNA.out.lsu_fasta
        .map { meta, fasta ->
            tuple(meta + [target: 'LSU'], fasta)
        }

    ch_five_eightS_fasta = DETECT_RNA.out.five_eightS_fasta
        .map { meta, fasta ->
            tuple(meta + [target: 'ITS'], fasta)
        }

    def all_fasta = ch_ssu_fasta.mix(ch_lsu_fasta,
                                     ch_five_eightS_fasta)

    mapseq_otu_dbs_in = mapseq_dbs_in.filter{ meta, _db -> meta.run_otu }
    MAPSEQ_OTU_KRONA(
        all_fasta,
        mapseq_otu_dbs_in
    )
    ch_versions = ch_versions.mix(MAPSEQ_OTU_KRONA.out.versions)

    // Infer amplified variable regions for SSU, extract reads for each amplified region if there are more than one //
    AMP_REGION_INFERENCE(
        DETECT_RNA.out.cmsearch_deoverlap_non_empty,
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
        DETECT_RNA.out.cmsearch_deoverlap_non_empty
    )
    ch_versions = ch_versions.mix(DADA2_SWF.out.versions)


    // ASV taxonomic assignments + generate Krona plots for each run+amp_region //
    mapseq_asv_dbs_in = mapseq_dbs_in.filter{ meta, _db -> meta.run_asv }
    MAPSEQ_ASV_KRONA(
        DADA2_SWF.out.dada2_out,
        AMP_REGION_INFERENCE.out.concat_var_regions,
        AMP_REGION_INFERENCE.out.extracted_var_path,
        mapseq_asv_dbs_in,
    )
    ch_versions = ch_versions.mix(MAPSEQ_ASV_KRONA.out.versions)

    // Count reads per ASV per var_region directly from DADA2 map files, without taxonomy //
    asv_read_counts_input = DADA2_SWF.out.dada2_out
        .map { meta, maps, _asvs, filt_reads ->
            [meta.subMap('id', 'single_end', 'var_regions_size'), meta.var_region, maps, filt_reads]
        }
        .transpose(by: 1)
        .join(
            AMP_REGION_INFERENCE.out.extracted_var_path
                .map { meta, var_region, filter_list ->
                    [meta.subMap('id', 'single_end', 'var_regions_size'), var_region, filter_list]
                },
            by: [0, 1]
        )
        .map { meta, var_region, maps, filt_reads, filter_list ->
            [meta + [var_region: var_region], maps, filt_reads, filter_list]
        }

    // Also count across all var regions combined for multi-region samples //
    concat_asv_read_counts_input = DADA2_SWF.out.dada2_out
        .map { meta, maps, _asvs, filt_reads ->
            [meta.subMap('id', 'single_end'), maps, filt_reads]
        }
        .join(
            AMP_REGION_INFERENCE.out.concat_var_regions
                .map { meta, concat_str, concat_vars ->
                    [meta.subMap('id', 'single_end'), concat_str, concat_vars]
                },
            by: [0]
        )
        .map { meta, maps, filt_reads, concat_str, concat_vars ->
            [meta + [var_region: concat_str], maps, filt_reads, concat_vars]
        }

    MAKE_ASV_COUNT_TABLES(asv_read_counts_input.mix(concat_asv_read_counts_input))
    ch_versions = ch_versions.mix(MAKE_ASV_COUNT_TABLES.out.versions)

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
            def json_content = new groovy.json.JsonBuilder(collected_json_maps).toPrettyString() 
            return json_content
        }
        .collectFile(
            name: "primer_validation_summary.json", 
            storeDir: "${params.outdir}", 
            newLine: true, 
            cache: false
        )
    

    
    /*****************************/
    /* ITS sanity check */
    /****************************/
    // Keep empty mask files in this channel as well: an empty mask must produce
    // empty filtered mseq files, rather than allowing the original annotations
    // to bypass this filtering step.
    rfam_masked_reads_mseq = MASK_FASTA_SWF.out.masked_out

    its_mseq_to_filter = MAPSEQ_OTU_KRONA.out.mseq
        .filter { meta, _mseq ->
            meta.db_label in ['ITSoneDB', 'UNITE']
        }
        .map { meta, mseq ->
            tuple(meta.id, meta, mseq)
        }
        .combine(
            rfam_masked_reads_mseq
                .map { meta, masked_reads ->
                    tuple(meta.id, masked_reads)
                },
            by: 0
        )
        .map { _id, meta, mseq, masked_reads ->
            tuple(meta, mseq, masked_reads)
        }

    FILTER_MAPSEQ_BY_MASK(its_mseq_to_filter)

    filtered_its_mseq = FILTER_MAPSEQ_BY_MASK.out.mseq

    its_mseq = filtered_its_mseq
    .map { meta, mseq ->
        // Use meta.id+target for the grouping key. This keeps one group per run+target,
        // even when ITS dbs are launched several times on the same run id.
        tuple(meta.id, meta.target, [(meta.db_label): mseq])
    }
    .groupTuple(by: [0, 1])
    .map { id, target, results_list ->
        def results = [:]

        results_list.each { result ->
            result.each { label, mseq ->
                // Make sure only non-empty mapseq files are included
                if (mseq.readLines().size() > 1) {
                    if (!results.containsKey(label))
                        results[label] = mseq
                }
            }
        }
        tuple([id: id, target: target], results)
    }
    // Exclude empty annotation files
    .filter { meta, results ->
        !results.isEmpty()
    }

    // Merge all mseq and rfam results for the same sample and target into a
    // single map keyed by db label (e.g. { 'ITSoneDB': mseq, 'UNITE': mseq,
    // 'Rfam_SSU_LSU_5_8S': masked_reads })
    read_assignments = its_mseq
    .map { meta, results ->
        tuple(meta.id, meta, results)
    }
    // Rfam masking is produced once per run, so combine it with every
    // target-specific ITS group rather than joining duplicate run ids.
    .combine(
        rfam_masked_reads_mseq
            .map { meta, masked_reads ->
                tuple(meta.id, masked_reads)
            },
        by: 0
    )
    .map { id, meta, results, masked_reads ->
        results['Rfam_SSU_LSU_5_8S'] = masked_reads
        [meta, results]
    }

    ITS_SANITY_CHECKER(read_assignments)

    // Keep one target-aware representation of every sanity-check result.
    // ITSoneDB and UNITE results for the same run+target are checked together,
    // while different targets remain separate.
    its_sanity_pass_boolean = ITS_SANITY_CHECKER.out.its_sanity_check_out
        .splitJson()
        .map { meta, test_results ->
            def pass = (
                test_results["tax_assignment_count_test"] &&
                test_results["mapping_proportion_test"] &&
                test_results["rank_proportion_test"]
            )
            tuple(meta.id, meta.target, pass)
        }

    // Derive run-level status from the target-aware sanity results. A run is
    // considered to have ITS taxonomies when at least one target passes.
    real_its_runs = its_sanity_pass_boolean
        .filter { _id, _target, passed -> passed }
        .map { id, _target, _passed -> [id: id] }
        .unique()

    // Identify target-specific ITS sanity checks that failed. The later run-level
    // status logic suppresses this failure when another target has taxonomy results.
    its_sanity_check_fails = its_sanity_pass_boolean
        .filter { _id, _target, passed -> !passed }
        .map { id, _target, _failed -> ["${id}", "failed"] }
        .unique()

    /*****************************/
    /* Publish OTU results */
    /****************************/

    // Join all MAPSEQ_OTU_KRONA outputs per sample+db combination
    otu_all_results = MAPSEQ_OTU_KRONA.out.mseq
        .join(MAPSEQ_OTU_KRONA.out.krona_input)
        .join(MAPSEQ_OTU_KRONA.out.biom_out)
        .join(MAPSEQ_OTU_KRONA.out.html)
        .filter { _meta, _mseq, krona_input, _biom_out, _html -> (krona_input.size() > 0) }

    // Branch into ITS and non-ITS databases
    otu_branched = otu_all_results
        .branch { meta, _mseq, _krona_input, _biom_out, _html ->
            its: meta.db_label in ["ITSoneDB", "UNITE"]
            non_its: meta.db_label !in ["ITSoneDB", "UNITE"]
        }

    // Retain only passing target instances, then merge the surviving MAPseq
    // files separately by database before regenerating their taxonomy outputs.
    merged_its_mseq = filtered_its_mseq
        .map { meta, mseq ->
            tuple(meta.id, meta.target, meta, mseq)
        }
        .combine(its_sanity_pass_boolean, by: [0, 1])
        .filter { _id, _target, _meta, mseq, passed -> passed && mseq.readLines().size() > 1 }
        .map { _id, _target, meta, mseq, _passed ->
            [meta.id, meta.db_id, meta.db_label, meta.target, mseq]
        }
        .groupTuple(by: [0, 1, 2])
        .map { id, db_id, db_label, targets, mseq_files ->
            tuple(
                [id: id, db_id: db_id, db_label: db_label, targets: targets],
                mseq_files
            )
        }

    MERGE_MAPSEQ(merged_its_mseq)

    its_db_inputs = mapseq_dbs_in
        .filter { meta, _db_files -> meta.db_label in ['ITSoneDB', 'UNITE'] }
        .map { meta, db_files ->
            def (_fasta, _tax, otu, _mscluster, _label) = db_files
            tuple(meta.id, otu, meta.db_label)
        }

    // MAPSEQ2BIOM needs the .OTU file from the input db
    // this step only involves ITS dbs as SSU+LSU mseq files stayed the same
    merged_its_mapseq2biom_input = MERGE_MAPSEQ.out.mseq
        .map { meta, mseq -> tuple(meta.db_id, meta, mseq) }
        .join(its_db_inputs)
        .multiMap { db_id, meta, mseq, otu, db_label ->
            mseq_ch: [meta, mseq]
            db_ch: [otu, db_label]
        }

    MAPSEQ2BIOM_MERGED(
        merged_its_mapseq2biom_input.mseq_ch,
        merged_its_mapseq2biom_input.db_ch,
    )
    ch_versions = ch_versions.mix(MAPSEQ2BIOM_MERGED.out.versions)

    KRONA_KTIMPORTTEXT_MERGED(MAPSEQ2BIOM_MERGED.out.krona_input)
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT_MERGED.out.versions)

    merged_its_results = MERGE_MAPSEQ.out.mseq
        .join(MAPSEQ2BIOM_MERGED.out.krona_input)
        .join(MAPSEQ2BIOM_MERGED.out.biom_out)
        .join(KRONA_KTIMPORTTEXT_MERGED.out.html)
        .map { meta, mseq, krona_input, biom_out, html ->
            [meta, mseq, krona_input, biom_out, html]
        }

    // Publish masked ITS reads once if at least one target instance passes.
    // List of unique runs passing one or more ITS sanity checks
    its_run_sanity_pass = its_sanity_pass_boolean
        .filter { _id, _target, passed -> passed }
        .map { id, _target, _passed -> [id, true] }
        .unique()

    passed_masked_out = MASK_FASTA_SWF.out.masked_out
        .filter { _meta, masked_fasta -> (masked_fasta.size() > 0) }
        .map { meta, masked_fasta ->
            [meta.id, meta, masked_fasta]
        }
        .combine(its_run_sanity_pass, by: 0)
        .filter { _id_meta, _meta, _masked_fasta, passed -> passed }
        .map { _id_meta, meta, masked_fasta, _passed ->
            tuple(meta, masked_fasta)
        }

    PUBLISH_MASKED_ITS(passed_masked_out)

    // Publish all non-ITS results + ITS results that pass sanity check
    publish_otu_input = otu_branched.non_its.mix(merged_its_results)
    PUBLISH_OTU_RESULTS(publish_otu_input)
    ch_versions = ch_versions.mix(PUBLISH_OTU_RESULTS.out.versions.first())

    /*****************************/
    /* MultiQC reports */
    /****************************/

    // Version collating //
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // Shared metadata is staged into every per-run and study report.
    ch_methods_description = Channel.value(
        WorkflowAmpliconpipeline.methodsDescriptionText(
            workflow,
            file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true),
            params
        )
    ).collectFile(name: 'methods_description_mqc.yaml')

    ch_multiqc_base_files = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml
        .mix(ch_methods_description)
        .collect()

    multiqc_input = CONCAT_PRIMER_CUTADAPT.out.cutadapt_json
        .map{ meta, json ->
            [meta.subMap('id', 'single_end'), json]
        }
        .join(
            READS_QC_MERGE.out.fastp_summary_json.map { meta, json ->
                [meta.subMap('id', 'single_end'), json]
            },
            remainder:true
        )
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
        .map { meta, final_inputs -> [meta.id, meta, final_inputs] }
        .join(
            ITS_SANITY_CHECKER.out.its_sanity_check_out_mqc
                .map { meta, mqc -> [meta.id, mqc] }
                .groupTuple(),
            remainder: true
        )
        .map { _id, meta, final_inputs, its_sanity_check_mqc ->
            def all_inputs = final_inputs instanceof List ? final_inputs.collect() : [final_inputs]
            if (its_sanity_check_mqc) {
                def sanity_check_files = its_sanity_check_mqc instanceof List ? its_sanity_check_mqc : [its_sanity_check_mqc]
                all_inputs.addAll(sanity_check_files)
            }
            [meta, all_inputs]
        }
    
    
    // MultiQC for individual runs //
    MULTIQC_RUN(
        multiqc_input,
        ch_multiqc_base_files,
        params.multiqc_config,
        [],
        [],
        [],
        [],
    )

    // generate aggregate summary of all its sanity check outputs
    study_its_sanity_check_path = ITS_SANITY_CHECKER.out.its_sanity_check_out_mqc
        .map { _meta, its_sanity_check -> its_sanity_check}
        .collectFile(name: "study_its_sanity_check_mqc.tsv", keepHeader: true, cache: false)

    // MultiQC for study !! assuming we do not have multiple studies in one samplesheet !! //
    multiqc_study = multiqc_input
        .flatten()
        .collect()
        .map { item -> item.findAll { it -> !(it instanceof Map) } }
        .map { dataList ->
            // have to remove the individual ITS sanity check outputs before including the study-wide file
            def its_files_to_remove = dataList.findAll { file -> file.name.contains("its_sanity_check_mqc.tsv") }
            dataList -= its_files_to_remove
        }
        .mix(study_its_sanity_check_path)
        .collect()
        .map { dataList -> [['id': 'study_multiqc_report'], dataList] }

    MULTIQC_STUDY(
        multiqc_study,
        ch_multiqc_base_files,
        params.multiqc_config,
        [],
        [],
        [],
        [],
    )

    /*****************************/
    /* End of execution reports */
    /****************************/


    // Runs can fail ITS sanity checking, but still have taxonomy results for various reasons
    // So we can't just automatically assign runs that fail ITS sanity checking as failed runs
    // And we definitely don't want runs that are labeled as succeeding and failing at the same time
    // The next bits of processing will handle this by grabbing all runs that have any kind of
    // taxonomic assignment results (SSU/LSU/ITS/ASV) and these will make up the subset of
    // passed runs. We will then filter out any runs that fail ITS sanity checking but are
    // in this subset of successful runs.

    // label runs that have ITS taxonomies (succeed at ITS sanity checking)
    passed_its_checks = real_its_runs
        .map{ meta ->
            [meta, "has_its_taxonomies"]
        }

    // label runs that have SSU/LSU taxonomies (non-ITS OTU results)
    has_ssu_lsu_taxonomies = otu_branched.non_its
        .map { meta, _mseq, _krona_input, _biom_out, _html ->
            [meta.subMap('id'), "has_ssu_lsu_taxonomies"]
        }
        .unique()

    // get status of runs that reach DADA2 but might fail for quality reasons
    dada2_stats_fail = DADA2_SWF.out.dada2_stats_fail.map { meta, stats_fail ->
        def key = meta.subMap('id')
        return [key, ["stats_fail": stats_fail]]
    }

    dada2_stats_fail
        .map{ meta, stats_fail ->
            if (stats_fail["stats_fail"] == "true"){
                [meta, "dada2_stats_fail"]
            }
            else{
                [meta, "has_dada2_results"]
            }
        }
        .set{ has_dada2_results }

    // groups all runs that have some taxonomy results
    all_taxonomy_results = has_dada2_results
        .mix(passed_its_checks, has_ssu_lsu_taxonomies)
        .groupTuple()
        .branch {
            _meta, results ->
                invalid: results.size() == 1 && results[0] == "dada2_stats_fail"
                valid: true
            // to remove cases where the only result a run has is "dada2_stats_fail" i.e. no OTU results either
        }

    // Extract passed runs, describe whether those passed runs also ASV results //
    // Rules are:
    //      if you have DADA2 results and SSU/LSU taxonomy results, you have `all_results`
    //      if you don't have DADA2 results but have ITS/SSU/LSU results
    //          if you have the dada2_stats_fail status, you have `dada2_stats_fail`
    //          if you don't have the dada2_stats_fail status, it means you have `no_asvs`
    final_passed_runs = all_taxonomy_results.valid
        .map { meta, results ->
            if ("has_dada2_results" in results &&
                ("has_ssu_lsu_taxonomies" in results || "has_its_taxonomies" in results)) {
                return "${meta.id},all_results"
            }
            else if ("has_dada2_results" in results) {
                return "${meta.id},asvs_only"
            }
            else if ("has_its_taxonomies" in results || "has_ssu_lsu_taxonomies" in results) {
                if ("dada2_stats_fail" in results) {
                    return "${meta.id},dada2_stats_fail"
                }
                else {
                    return "${meta.id},no_asvs"
                }
            }
            error("Unexpected. meta: ${meta}, results: ${results}")
        }

    // Save all passed runs to file //
    final_passed_runs
        .collectFile(name: "qc_passed_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)

    runs_with_taxonomies = all_taxonomy_results.valid
        .map{ meta, _results -> ["${meta.id}", "passed"] }

    // Extract runs that failed SeqFu check //
    seqfu_fails = READS_QC.out.seqfu_check
        .splitCsv(sep: "\t", elem: 1)
        .filter { _meta, seqfu_res ->
            seqfu_res[0] != "OK"
        }
        .map { meta, __ -> "${meta.id},seqfu_fail" }

    // Extract runs that failed Suffix Header check //
    sfxhd_fails = READS_QC.out.suffix_header_check
        .filter { _meta, sfxhd_res ->
            sfxhd_res.countLines() != 0
        }
        .map { meta, __ -> "${meta.id},sfxhd_fail" }

    // Extract runs that failed Library Strategy check //
    libstrat_fails = READS_QC_MERGE.out.amplicon_check
        .filter { _meta, strategy ->
            strategy != "AMPLICON"
        }
        .map { meta, __ -> "${meta.id},libstrat_fail" }

    // Extract runs that had zero reads after fastp //
    no_reads_fails = extended_reads_qc.qc_empty
        .map { meta, __ -> "${meta.id},no_reads" }
        .mix(
            DETECT_RNA.out.cmsearch_deoverlap_empty
                .map { meta, __ -> "${meta.id},no_reads" }
        )

    // Extract runs that have no OTU results and a "dada2_stats_fail" status for ASVs
    no_otus_dada2_stats_fails = all_taxonomy_results.invalid
        .map { meta, __ -> "${meta.id}, no_otus_dada2_stats_fail"}

    // filter out runs that fail ITS sanity checking but have other taxonomy results
    failed_its_runs = its_sanity_check_fails
        .mix(runs_with_taxonomies)
        .groupTuple()
        .filter{ _run, result ->
            result == ["failed"]
        }
        .map { run, _result ->
            "${run},its_sanity_check_fail"
        }

    // Save all failed runs to file //
    seqfu_fails
        .concat(sfxhd_fails, libstrat_fails, no_reads_fails, failed_its_runs, no_otus_dada2_stats_fails)
        .collectFile(name: "qc_failed_runs.csv", storeDir: "${params.outdir}", newLine: true, cache: false)

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
                    "identification_strategy": name.contains("_auto") ? "auto" : "std",
                ]
                json_map["primers"] << new_primer
            }

            json_map
        }
        .collect()
        .map { collected_json_maps -> 
            def json_content = new groovy.json.JsonBuilder(collected_json_maps).toPrettyString() 
            return json_content
        }
        .collectFile(name: "primer_validation_summary.json", storeDir: "${params.outdir}", newLine: true, cache: false)
}
