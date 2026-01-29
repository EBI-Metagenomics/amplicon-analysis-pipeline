
include { MAPSEQ                } from '../../modules/ebi-metagenomics/mapseq/main'
include { MAPSEQ2ASVTABLE       } from '../../modules/local/mapseq2asvtable/main.nf'
include { MAKE_ASV_COUNT_TABLES } from '../../modules/local/make_asv_count_tables/main.nf'
include { KRONA_KTIMPORTTEXT    } from '../../modules/ebi-metagenomics/krona/ktimporttext/main'
include { EXTRACT_ASVS_LEFT     } from '../../modules/local/extract_asvs_left/main.nf'

workflow MAPSEQ_ASV_KRONA {
    
    take:
        dada2_output
        concat_var_regions
        extracted_var_path
        dbs_in

    main:

        ch_versions = channel.empty()

        mapseq_in = dada2_output
            .map{ meta, maps, asv_seqs, filt_reads -> [meta, [maps, asv_seqs, filt_reads]] }
            .combine(dbs_in)
            .map { asv_meta, asv_files, db_meta, db_files ->
                def meta = asv_meta + ['db_id': db_meta.id, 'db_label': db_meta.label, 'dada2_label': db_meta.dada2_label]
                def (fasta, tax, _otu, mscluster, label) = db_files
                def (_maps, asv_seqs, _filt_reads) = asv_files
                return [meta.id, meta, asv_seqs, fasta, tax, mscluster, label]
            }
            .groupTuple()
            .map { _meta_id, vs -> [vs.size(), vs] }
            .flatten()
            .map{ n, meta, asv_seqs, fasta, tax, mscluster, label -> 
                  [groupKey(meta, size: n), asv_seqs, fasta, tax, mscluster, label] }

        MAPSEQ(mapseq_in)
        ch_versions = ch_versions.mix(MAPSEQ.out.versions.first())

        mapseq2biom_in = MAPSEQ.out.mseq
            .map { meta, mapseq_out -> [meta, mapseq_out, meta.dada2_label] }
        MAPSEQ2ASVTABLE(mapseq2biom_in)
        ch_versions = ch_versions.mix(MAPSEQ2ASVTABLE.out.versions.first())

        // Transpose by var region in case any samples have more than one
        split_mapseq2asvtable = MAPSEQ2ASVTABLE.out.asvtaxtable
            .groupTuple()
            .map { meta, asvtaxtables -> 
                   [meta.subMap('id', 'single_end', 'var_regions_size'), meta['var_region'], meta.dada2_label, asvtaxtables] }
            .transpose(by: 1)
        split_mapseq2asvtable.view{ it -> "split_mapseq2asvtable - ${it}"}


        // Transpose by var region in case any samples have more than one. Also reorder the inputs slightly
        split_input_ = dada2_output
            .map { meta, maps, _asv_seqs, filt_reads -> 
                   [meta.subMap('id', 'single_end', 'var_regions_size'), meta['var_region'], maps, filt_reads] }
            .transpose(by: 1)
            .join(extracted_var_path, by: [0, 1])
            .join(split_mapseq2asvtable, by: [0, 1])  // here's my trouble, could join and then flatten, or use a cross
        split_input_.view{ it -> "split_input_ - ${it}"}
        
        split_input = split_input_
            .transpose(by: -1)
            .map { _submeta, var_region, meta, maps, filt_reads, extracted_var, dada2_label, asvtaxtable ->
                   [meta, var_region, dada2_label, maps, asvtaxtable, filt_reads, extracted_var] }
        split_input.view{ it -> "split_input - ${it}"}

        // Make a channel containing the concatenated var region for any sample that has more than one var region
        multi_region_concats = split_input
            .map { meta, var_region, dada2_label, maps, asvtaxtable, filt_reads, extracted_var ->
                   [meta.subMap('id', 'single_end'), var_region, dada2_label, meta['var_regions_size'], maps, asvtaxtable, filt_reads, extracted_var] }
            .join(concat_var_regions, by: 0)
            .map { meta, _var_region, dada2_label, var_regions_size, maps, asvtaxtable, filt_reads, _extracted_vars, concat_str, concat_vars ->
                   [meta + ['var_regions_size':var_regions_size], concat_str, dada2_label, maps, asvtaxtable, filt_reads, concat_vars] }

        // Add in the concatenated var region channel to the rest of the input
        final_asv_count_table_input = split_input
            .mix(multi_region_concats)
            .map { meta, var_region, dada2_label, maps, asvtaxtable, filt_reads, extracted_var ->
                   [meta + ['var_region': var_region], maps, asvtaxtable, filt_reads, extracted_var, dada2_label] }
        
        MAKE_ASV_COUNT_TABLES(final_asv_count_table_input)
        ch_versions = ch_versions.mix(MAKE_ASV_COUNT_TABLES.out.versions.first())

        KRONA_KTIMPORTTEXT(
            MAKE_ASV_COUNT_TABLES.out.asv_krona_counts,
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions.first())

    emit:
        asv_krona_counts = MAKE_ASV_COUNT_TABLES.out.asv_krona_counts
        asv_read_counts = MAKE_ASV_COUNT_TABLES.out.asv_read_counts
        asvtaxtable = MAPSEQ2ASVTABLE.out.asvtaxtable
        krona_out = KRONA_KTIMPORTTEXT.out.html
        versions = ch_versions
}
