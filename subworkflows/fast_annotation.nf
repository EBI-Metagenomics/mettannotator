include { AMRFINDER_PLUS_GET_SPECIES_NAME            } from '../modules/local/amrfinder_plus'
include { AMRFINDER_PLUS; AMRFINDER_PLUS_TO_GFF      } from '../modules/local/amrfinder_plus'
include { AMRFINDER_PLUS_TSV_POSTPROCESSING          } from '../modules/local/amrfinder_plus'
include { DEFENSE_FINDER                             } from '../modules/local/defense_finder'
include { CRISPRCAS_FINDER                           } from '../modules/local/crisprcasfinder'
include { EGGNOG_MAPPER as EGGNOG_MAPPER_ORTHOLOGS   } from '../modules/local/eggnog'
include { EGGNOG_MAPPER as EGGNOG_MAPPER_ANNOTATIONS } from '../modules/local/eggnog'
include { DETECT_TRNA                                } from '../modules/local/detect_trna'
include { DETECT_NCRNA                               } from '../modules/local/detect_ncrna'
include { ANTISMASH                                  } from '../modules/local/antismash'
include { ANTISMASH_TO_GFF                           } from '../modules/local/antismash'
include { ANTISMASH_SUMMARY                          } from '../modules/local/antismash'
include { DBCAN                                      } from '../modules/local/dbcan'
include { PSEUDOFINDER                               } from '../modules/local/pseudofinder'
include { PSEUDOFINDER_POSTPROCESSING                } from '../modules/local/pseudofinder'
include { SELECT_HYPOTHETICAL_PROTEINS               } from '../modules/local/select_hypothetical_proteins'
include { DIAMOND_BLASTP as DIAMOND_UNIREF90         } from '../modules/local/diamond_blastp'
include { DIAMOND_BLASTP as DIAMOND_UNIREF50         } from '../modules/local/diamond_blastp'
include { GECCO_RUN                                  } from '../modules/nf-core/gecco/run/main'
include { QUAST                                      } from '../modules/nf-core/quast/main'

workflow FAST_ANNOTATION {

    take:
    faa              // channel: [meta, faa]
    gff              // channel: [meta, gff]
    gbk              // channel: [meta, gbk]
    fna              // channel: [meta, fna]
    compliant_gbk    // channel: [meta, compliant_gbk]
    compliant_gff    // channel: [meta, compliant_gff]
    assemblies       // channel: [meta, assembly]
    detected_kingdom // channel: [meta, kingdom]
    amrfinder_plus_db
    antismash_db
    defense_finder_db
    dbcan_db
    eggnog_db
    rfam_ncrna_models
    pseudofinder_db
    uniref90_db      // channel: [path, version]
    uniref50_db      // channel: [path, version]

    main:
    ch_versions = Channel.empty()

    assemblies_and_gff = assemblies.join(gff)
    QUAST(
        assemblies_and_gff.map { meta, assembly, _gff     -> [meta, assembly]  },
        assemblies_and_gff.map { meta, _assembly, gff_file -> [meta, gff_file] }
    )
    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    PSEUDOFINDER(compliant_gbk, pseudofinder_db)
    ch_versions = ch_versions.mix(PSEUDOFINDER.out.versions.first())

    PSEUDOFINDER_POSTPROCESSING(
        gff.join(compliant_gff).join(PSEUDOFINDER.out.pseudofinder_gff)
    )
    ch_versions = ch_versions.mix(PSEUDOFINDER_POSTPROCESSING.out.versions.first())

    CRISPRCAS_FINDER(assemblies)
    ch_versions = ch_versions.mix(CRISPRCAS_FINDER.out.versions.first())

    proteins_for_emapper_orth = faa.map { meta, faa_file -> tuple(meta, file(faa_file), file("NO_FILE")) }
    EGGNOG_MAPPER_ORTHOLOGS(proteins_for_emapper_orth, Channel.value("mapper"), eggnog_db)
    ch_versions = ch_versions.mix(EGGNOG_MAPPER_ORTHOLOGS.out.versions.first())

    orthologs_for_annotations = assemblies.join(EGGNOG_MAPPER_ORTHOLOGS.out.orthologs).map {
        meta, _assembly, orthologs -> tuple(meta, file("NO_FILE"), file(orthologs))
    }
    EGGNOG_MAPPER_ANNOTATIONS(orthologs_for_annotations, Channel.value("annotations"), eggnog_db)
    ch_versions = ch_versions.mix(EGGNOG_MAPPER_ANNOTATIONS.out.versions.first())

    AMRFINDER_PLUS_GET_SPECIES_NAME(assemblies, amrfinder_plus_db)

    AMRFINDER_PLUS(
        assemblies.join(faa).join(gff).join(AMRFINDER_PLUS_GET_SPECIES_NAME.out.detected_organism),
        amrfinder_plus_db
    )
    ch_versions = ch_versions.mix(AMRFINDER_PLUS.out.versions.first())

    AMRFINDER_PLUS_TSV_POSTPROCESSING(AMRFINDER_PLUS.out.amrfinder_tsv)
    ch_versions = ch_versions.mix(AMRFINDER_PLUS_TSV_POSTPROCESSING.out.versions.first())

    AMRFINDER_PLUS_TO_GFF(AMRFINDER_PLUS_TSV_POSTPROCESSING.out.amrfinder_tsv)
    ch_versions = ch_versions.mix(AMRFINDER_PLUS_TO_GFF.out.versions.first())

    DEFENSE_FINDER(faa.join(gff), defense_finder_db)
    ch_versions = ch_versions.mix(DEFENSE_FINDER.out.versions.first())

    DETECT_TRNA(fna.join(detected_kingdom))
    ch_versions = ch_versions.mix(DETECT_TRNA.out.versions.first())

    DETECT_NCRNA(fna, rfam_ncrna_models)
    ch_versions = ch_versions.mix(DETECT_NCRNA.out.versions.first())

    GECCO_RUN(gbk.map { meta, gbk_file -> [meta, gbk_file, []] }, [])
    ch_versions = ch_versions.mix(GECCO_RUN.out.versions.first())

    ANTISMASH(gbk, antismash_db)
    ch_versions = ch_versions.mix(ANTISMASH.out.versions.first())

    ANTISMASH_TO_GFF(ANTISMASH.out.json)
    ch_versions = ch_versions.mix(ANTISMASH_TO_GFF.out.versions.first())

    ANTISMASH_SUMMARY(ANTISMASH_TO_GFF.out.gff)
    ch_versions = ch_versions.mix(ANTISMASH_SUMMARY.out.versions.first())

    DBCAN(faa.join(gff), dbcan_db)
    ch_versions = ch_versions.mix(DBCAN.out.versions.first())

    SELECT_HYPOTHETICAL_PROTEINS(faa)

    DIAMOND_UNIREF90(SELECT_HYPOTHETICAL_PROTEINS.out.hypothetical_faa, uniref90_db, "uniref90")
    ch_versions = ch_versions.mix(DIAMOND_UNIREF90.out.versions.first())

    DIAMOND_UNIREF50(SELECT_HYPOTHETICAL_PROTEINS.out.hypothetical_faa, uniref50_db, "uniref50")
    ch_versions = ch_versions.mix(DIAMOND_UNIREF50.out.versions.first())

    emit:
    eggnog_annotations = EGGNOG_MAPPER_ANNOTATIONS.out.annotations
    ncrna_tblout       = DETECT_NCRNA.out.ncrna_tblout
    trna_gff           = DETECT_TRNA.out.trna_gff
    crisprcas_hq_gff   = CRISPRCAS_FINDER.out.hq_gff
    amrfinder_tsv      = AMRFINDER_PLUS.out.amrfinder_tsv
    antismash_gff      = ANTISMASH_TO_GFF.out.gff
    gecco_gff          = GECCO_RUN.out.gff
    dbcan_gff          = DBCAN.out.dbcan_gff
    defense_gff        = DEFENSE_FINDER.out.gff
    pseudofinder_gff   = PSEUDOFINDER_POSTPROCESSING.out.pseudofinder_processed_gff
    quast_results      = QUAST.out.results
    uniref90_tsv       = DIAMOND_UNIREF90.out.hits_tsv
    uniref50_tsv       = DIAMOND_UNIREF50.out.hits_tsv
    versions           = ch_versions
}
