
def dada2_input_preparation_function( concat_input, reads_qc, cutadapt_channel ) {

    def dada2_input = concat_input
        // Stage 1: Collapse per-region primer results into one tuple per sample.
        // Uses groupKey sized by var_regions_size so groupTuple emits early.
        // See: https://training.nextflow.io/advanced/grouping/#grouping-using-submap
        .map { meta, _std_primer, _auto_primer ->
            def key = groupKey(meta.subMap('id', 'single_end'), meta['var_regions_size'])
            [ key, meta['var_region'], meta['var_regions_size'] ]
        }
        .groupTuple(by: 0)

        // Stage 2: Pair each sample with both its fastp reads (via join) and
        // cutadapt reads (via mix), then group them together using groupKey of
        // size 2 so we get both read sources per sample.
        // Note: reads[0] = fastp, reads[1] = cutadapt after groupTuple.
        .join(
            reads_qc.map{ meta, reads -> [meta.subMap('id', 'single_end'), reads] }, 
            by: 0
        )
        .mix(cutadapt_channel)
        .map { meta, var_region, var_regions_size, reads ->
            [ groupKey(meta.subMap('id', 'single_end'), 2), var_region, var_regions_size, reads ]
        }
        .groupTuple(by: 0)

        // Stage 3: Reconstruct meta and select reads — prefer cutadapt output
        // (primer-trimmed) if non-empty, otherwise fall back to fastp output.
        .map { meta, var_region, var_regions_size, reads ->
            def final_meta = meta + [
                'var_region':      var_region.unique()[0],
                'var_regions_size': var_regions_size[0][0]
            ]
            def fastp_reads = reads[0]
            def cutadapt_reads = reads[1]
            // For SE reads, .size() returns file size directly;
            // for PE reads (a list), check the first file's size.
            def cutadapt_read_size = meta.single_end
                ? cutadapt_reads.size()
                : cutadapt_reads[0].size()
            def final_reads = cutadapt_read_size > 0 ? cutadapt_reads : fastp_reads
            [ final_meta, final_reads ]
        }

    return dada2_input
}
