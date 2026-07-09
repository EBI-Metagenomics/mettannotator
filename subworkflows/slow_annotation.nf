include { ADD_TAXID_TO_PROTEIN_FASTA } from '../modules/local/add_taxid'
include { INTERPROSCAN               } from '../modules/local/interproscan'
include { UNIFIRE                    } from '../modules/local/unifire'
include { SANNTIS                    } from '../modules/local/sanntis'

workflow SLOW_ANNOTATION {

    take:
    faa             // channel: [meta, faa]
    gbk             // channel: [meta, gbk]
    interproscan_db // channel: path

    main:
    ch_versions = Channel.empty()

    ADD_TAXID_TO_PROTEIN_FASTA(faa)
    ch_versions = ch_versions.mix(ADD_TAXID_TO_PROTEIN_FASTA.out.versions.first())

    INTERPROSCAN(ADD_TAXID_TO_PROTEIN_FASTA.out.annotations_faa_with_taxid, interproscan_db)
    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions.first())

    UNIFIRE(INTERPROSCAN.out.ips_xml)
    ch_versions = ch_versions.mix(UNIFIRE.out.versions.first())

    SANNTIS(INTERPROSCAN.out.ips_annotations.join(gbk))
    ch_versions = ch_versions.mix(SANNTIS.out.versions.first())

    emit:
    ips_annotations = INTERPROSCAN.out.ips_annotations
    sanntis_gff     = SANNTIS.out.sanntis_gff
    arba            = UNIFIRE.out.arba
    unirule         = UNIFIRE.out.unirule
    pirsr           = UNIFIRE.out.pirsr
    versions        = ch_versions
}
