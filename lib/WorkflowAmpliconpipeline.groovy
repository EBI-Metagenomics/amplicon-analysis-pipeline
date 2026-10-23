//
// This file holds several functions specific to the workflow/ampliconpipeline.nf in the nf-core/ampliconpipeline pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowAmpliconpipeline {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        genomeExistsError(params, log)


        if (!params.fasta) {
            Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        def citation_text = [
            "The pipeline supports read quality control with fastp (Chen et al. 2018),",
            "primer trimming with Cutadapt (Martin 2011),",
            "RNA homology searches with Infernal (Nawrocki &amp; Eddy 2013),",
            "amplicon sequence variant inference with DADA2 (Callahan et al. 2016),",
            "taxonomic classification with MAPseq (Matias Rodrigues et al. 2017),",
            "taxonomic visualisation with Krona (Ondov et al. 2011),",
            "and report aggregation with MultiQC (Ewels et al. 2016).",
            "Supporting steps use SeqFu for FASTQ validation (Telatin et al. 2021),",
            "SeqKit for sequence processing (Shen et al. 2024),",
            "BEDTools for sequence masking (Quinlan &amp; Hall 2010),",
            "PIMENTO for primer inference (Atallah et al. 2026),",
            "and BBMap for optional read standardisation (Bushnell 2014).",
            "Sequence extraction and result processing also use seqtk, Easel (Eddy 2011),",
            "cmsearch_tblout_deoverlap, and mgnify-pipelines-toolkit."
        ].join(' ')

        return citation_text
    }

    public static String toolBibliographyText(params) {

        def reference_text = [
            '<li>Chen, S., Zhou, Y., Chen, Y., &amp; Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34, i884–i890. doi: <a href="https://doi.org/10.1093/bioinformatics/bty560">10.1093/bioinformatics/bty560</a></li>',
            '<li>Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17, 10–12. doi: <a href="https://doi.org/10.14806/ej.17.1.200">10.14806/ej.17.1.200</a></li>',
            '<li>Nawrocki, E. P., &amp; Eddy, S. R. (2013). Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics, 29, 2933–2935. doi: <a href="https://doi.org/10.1093/bioinformatics/btt509">10.1093/bioinformatics/btt509</a></li>',
            '<li>Callahan, B. J., et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. Nature Methods, 13, 581–583. doi: <a href="https://doi.org/10.1038/nmeth.3869">10.1038/nmeth.3869</a></li>',
            '<li>Matias Rodrigues, J. F., Schmidt, T. S. B., Tackmann, J., &amp; von Mering, C. (2017). MAPseq: highly efficient k-mer search with confidence estimates, for rRNA sequence analysis. Bioinformatics, 33, 3808–3810. doi: <a href="https://doi.org/10.1093/bioinformatics/btx517">10.1093/bioinformatics/btx517</a></li>',
            '<li>Ondov, B. D., Bergman, N. H., &amp; Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 12, 385. doi: <a href="https://doi.org/10.1186/1471-2105-12-385">10.1186/1471-2105-12-385</a></li>',
            '<li>Ewels, P., Magnusson, M., Lundin, S., &amp; Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32, 3047–3048. doi: <a href="https://doi.org/10.1093/bioinformatics/btw354">10.1093/bioinformatics/btw354</a></li>',
            '<li>Telatin, A., Fariselli, P., &amp; Birolo, G. (2021). SeqFu: A Suite of Utilities for the Robust and Reproducible Manipulation of Sequence Files. Bioengineering, 8(5), 59. <a href="https://doi.org/10.3390/bioengineering8050059">SeqFu</a>.</li>',
            '<li>Shen, W., Sipos, B., &amp; Zhao, L. (2024). SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing. iMeta, e191. <a href="https://doi.org/10.1002/imt2.191">SeqKit</a>.</li>',
            '<li>Li, H. seqtk: toolkit for processing sequences in FASTA/Q formats. Software project. <a href="https://github.com/lh3/seqtk">seqtk</a>.</li>',
            '<li>Quinlan, A. R., &amp; Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841–842. <a href="https://doi.org/10.1093/bioinformatics/btq033">BEDTools</a>.</li>',
            '<li>Eddy, S. R. (2011). Accelerated Profile HMM Searches. PLoS Computational Biology, 7(10), e1002195. doi: <a href="https://doi.org/10.1371/journal.pcbi.1002195">10.1371/journal.pcbi.1002195</a> (Easel / HMMER).</li>',
            '<li>Nawrocki, E. cmsearch_tblout_deoverlap: remove lower-scoring overlaps from Infernal cmsearch results. Software project. <a href="https://github.com/nawrockie/cmsearch_tblout_deoverlap">cmsearch_tblout_deoverlap</a>.</li>',
            '<li>Atallah, C., Richardson, L., Beracochea, M., &amp; Finn, R. D. (2026). PIMENTO: a primer inference toolkit to facilitate large-scale calling of amplicon sequence variants. GigaScience, giag083. <a href="https://doi.org/10.1093/gigascience/giag083">PIMENTO</a>.</li>',
            '<li>EMBL-EBI Microbiome Informatics Team. mgnify-pipelines-toolkit: scripts and tools for MGnify pipelines. Software project. <a href="https://github.com/EBI-Metagenomics/mgnify-pipelines-toolkit">mgnify-pipelines-toolkit</a>.</li>',
            '<li>Bushnell, B. (2014). BBMap: A Fast, Accurate, Splice-Aware Aligner. Lawrence Berkeley National Laboratory, LBNL-7065E. <a href="https://escholarship.org/uc/item/1h3515gn">BBMap</a>.</li>'
        ].join('\n  ')

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Core capabilities are described here; the software versions section records tools used in this run.
        meta["tool_citations"] = toolCitationText(params)
        meta["tool_bibliography"] = toolBibliographyText(params)

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html.toString()
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }
}
