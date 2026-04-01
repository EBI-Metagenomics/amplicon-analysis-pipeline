
include { MAPSEQ             } from '../../../modules/ebi-metagenomics/mapseq/main'
include { MAPSEQ2BIOM        } from '../../../modules/ebi-metagenomics/mapseq2biom/main'
include { KRONA_KTIMPORTTEXT } from '../../../modules/ebi-metagenomics/krona/ktimporttext/main'

workflow MAPSEQ_OTU_KRONA {

    take:
    ch_fasta    // channel: [ val(meta), [ fasta ] ]
    ch_dbs      // channel: [ val(meta), [ path(fasta), path(tax), path(otu), path(mscluster), val(label) ] ]

    main:

    ch_versions = channel.empty()

    input = ch_fasta
        .combine(ch_dbs)
        .filter { reads_meta, _reads, db_meta, _db_files ->
            reads_meta.db_id == null || reads_meta.db_id == db_meta.id
        }
        .map { reads_meta, reads, db_meta, db_files ->
            def meta = reads_meta + ['db_id': db_meta.id, 'db_label': db_meta.db_label]
            def (fasta, tax, otu, mscluster, label) = db_files
            return [meta, reads, fasta, tax, otu, mscluster, label]
        }

    mapseq_in = input
        .multiMap { meta, reads, fasta, tax, _otu, mscluster, _label ->
            reads_ch: [meta, reads]
            db_ch: [fasta, tax, mscluster]
        }

    MAPSEQ(mapseq_in.reads_ch, mapseq_in.db_ch)
    ch_versions = ch_versions.mix(MAPSEQ.out.versions.first())

    mapseq2biom_in = MAPSEQ.out.mseq
        .join(input)
        .multiMap { meta, mapseq_out, _reads, _fasta, _tax, otu, _mscluster, label ->
            mseq_ch: [meta, mapseq_out]
            db_ch: [otu, label]
        }

    MAPSEQ2BIOM(mapseq2biom_in.mseq_ch, mapseq2biom_in.db_ch)
    ch_versions = ch_versions.mix(MAPSEQ2BIOM.out.versions.first())

    KRONA_KTIMPORTTEXT(MAPSEQ2BIOM.out.krona_input)
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions.first())

    emit:
    mseq             = MAPSEQ.out.mseq                      // channel: [ val(meta), [ mseq ] ]
    krona_input      = MAPSEQ2BIOM.out.krona_input         // channel: [ val(meta), [ txt ] ]
    biom_out         = MAPSEQ2BIOM.out.biom_out            // channel: [ val(meta), [ tsv ] ]
    biom_notaxid_out = MAPSEQ2BIOM.out.biom_notaxid_out    // channel: [ val(meta), [ tsv ] ]
    html             = KRONA_KTIMPORTTEXT.out.html         // channel: [ val(meta), [ html ] ]
    versions         = ch_versions                         // channel: [ versions.yml ]
}
