include { validateParameters; paramsHelp; paramsSummaryLog; } from 'plugin/nf-schema'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run main.nf --input samplesheet.csv --outdir ./outdir")
   exit 0
}
// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

include { AMPLICON_PIPELINE } from './workflows/pipeline.nf'

workflow {
    AMPLICON_PIPELINE()
}