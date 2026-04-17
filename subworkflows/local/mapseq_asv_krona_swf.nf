
include { MAPSEQ                } from '../../modules/ebi-metagenomics/mapseq/main'
include { MAPSEQ2ASVTABLE       } from '../../modules/local/mapseq2asvtable/main.nf'
include { MAKE_ASV_COUNT_TABLES } from '../../modules/local/make_asv_count_tables/main.nf'
include { KRONA_KTIMPORTTEXT    } from '../../modules/ebi-metagenomics/krona/ktimporttext/main'

workflow MAPSEQ_ASV_KRONA {
    
    take:
        dada2_output
        concat_var_regions
        extracted_var_path
        dbs_in

    main:

        ch_versions = channel.empty()

        mapseq_in = dada2_output
            .map{ meta, maps, asv_seqs, _filt_reads -> [meta, [maps, asv_seqs]] }
            .combine(dbs_in)
            .map { asv_meta, asv_files, db_meta, db_files ->
                def meta = asv_meta + ['db_id': db_meta.id, 'asv_label': db_meta.asv_label]
                def (fasta, tax, _otu, mscluster, _label) = db_files
                def (_maps, asv_seqs) = asv_files
                return [meta.subMap('id', 'single_end', 'var_region', 'var_regions_size'), [meta, asv_seqs, fasta, tax, mscluster, meta.asv_label]]
            }
            // Group (groupTuple()) to count number of items in each group, to be used in a future groupKey
            .groupTuple()
            .map { meta_k, vs -> [meta_k, vs.size(), vs] }
            // then re-expand (transpose(by: 2)) with each item now carrying n in its meta.
            .transpose(by: 2)
            .map{ _meta_k, n, v ->
                def (meta, asv_seqs, fasta, tax, mscluster, label) = v
                return [meta + ['n': n], asv_seqs, fasta, tax, mscluster, label]
            }

        mapseq_in
            .multiMap { meta, asv_seqs, fasta, tax, mscluster, _label ->
                reads_ch: [meta, asv_seqs]
                db_ch: [fasta, tax, mscluster]
            }
            .set { mapseq_split }

        MAPSEQ(mapseq_split.reads_ch, mapseq_split.db_ch)
        ch_versions = ch_versions.mix(MAPSEQ.out.versions.first())

        mapseq2asvtable_in = MAPSEQ.out.mseq
            .map { meta, mapseq_out -> [meta, mapseq_out, meta.asv_label] }
        MAPSEQ2ASVTABLE(mapseq2asvtable_in)
        ch_versions = ch_versions.mix(MAPSEQ2ASVTABLE.out.versions.first())

        // Transpose by var region in case any samples have more than one
        split_mapseq2asvtable = MAPSEQ2ASVTABLE.out.asvtaxtable
            .map{ meta, asvtaxtable -> 
                  [ groupKey(meta.subMap('id', 'single_end', 'var_region', 'var_regions_size'), meta.n), [meta.asv_label, asvtaxtable] ] }
            .groupTuple()
            .map { meta, asvtaxtables -> 
                   [meta.subMap('id', 'single_end', 'var_regions_size'), meta['var_region'], asvtaxtables] }
            .transpose(by: 1)

        // Transpose by var region in case any samples have more than one. Also reorder the inputs slightly
        split_input = dada2_output
            .map { meta, maps, _asv_seqs, filt_reads -> 
                   [meta.subMap('id', 'single_end', 'var_regions_size'), meta['var_region'], maps, filt_reads] }
            .transpose(by: 1)
            .join(extracted_var_path, by: [0, 1])
            .join(split_mapseq2asvtable, by: [0, 1])
            .transpose(by: 5)
            .map { meta, var_region, maps, filt_reads, extracted_var, asvlabel_asvtaxtable ->
                def (asv_label, asvtaxtable) = asvlabel_asvtaxtable.flatten()
                return [meta, var_region, asv_label, maps, asvtaxtable, filt_reads, extracted_var]
            }

        // Make a channel containing the concatenated var region for any sample that has more than one var region
        multi_region_concats = split_input
            .map { meta, var_region, asv_label, maps, asvtaxtable, filt_reads, extracted_var ->
                   [meta.subMap('id', 'single_end'), var_region, asv_label, meta['var_regions_size'], maps, asvtaxtable, filt_reads, extracted_var] }
            .combine(concat_var_regions, by: 0)
            .map { meta, _var_region, asv_label, var_regions_size, maps, asvtaxtable, filt_reads, _extracted_vars, concat_str, concat_vars ->
                   [meta + ['var_regions_size':var_regions_size], concat_str, asv_label, maps, asvtaxtable, filt_reads, concat_vars] }

        // Add in the concatenated var region channel to the rest of the input
        final_asv_count_table_input = split_input
            .mix(multi_region_concats)
            .map { meta, var_region, asv_label, maps, asvtaxtable, filt_reads, extracted_var ->
                   [meta + ['var_region': var_region, 'asv_label': asv_label], maps, asvtaxtable, filt_reads, extracted_var, asv_label] }

        MAKE_ASV_COUNT_TABLES(final_asv_count_table_input)
        ch_versions = ch_versions.mix(MAKE_ASV_COUNT_TABLES.out.versions.first())

        KRONA_KTIMPORTTEXT(MAKE_ASV_COUNT_TABLES.out.asv_krona_counts)
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions.first())

    emit:
        asv_krona_counts = MAKE_ASV_COUNT_TABLES.out.asv_krona_counts
        asv_read_counts = MAKE_ASV_COUNT_TABLES.out.asv_read_counts
        asvtaxtable = MAPSEQ2ASVTABLE.out.asvtaxtable
        krona_out = KRONA_KTIMPORTTEXT.out.html
        versions = ch_versions
}
