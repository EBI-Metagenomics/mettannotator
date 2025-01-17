/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { LOOKUP_KINGDOM                             } from '../modules/local/lookup_kingdom'
include { PROKKA as PROKKA_STANDARD                  } from '../modules/local/prokka'
include { PROKKA as PROKKA_COMPLIANT                 } from '../modules/local/prokka'
include { AMRFINDER_PLUS; AMRFINDER_PLUS_TO_GFF      } from '../modules/local/amrfinder_plus'
include { DEFENSE_FINDER                             } from '../modules/local/defense_finder'
include { CRISPRCAS_FINDER                           } from '../modules/local/crisprcasfinder'
include { EGGNOG_MAPPER as EGGNOG_MAPPER_ORTHOLOGS   } from '../modules/local/eggnog'
include { EGGNOG_MAPPER as EGGNOG_MAPPER_ANNOTATIONS } from '../modules/local/eggnog'
include { INTERPROSCAN                               } from '../modules/local/interproscan'
include { DETECT_TRNA                                } from '../modules/local/detect_trna'
include { DETECT_NCRNA                               } from '../modules/local/detect_ncrna'
include { SANNTIS                                    } from '../modules/local/sanntis'
include { UNIFIRE                                    } from '../modules/local/unifire'
include { ANNOTATE_GFF                               } from '../modules/local/annotate_gff'
include { ANTISMASH                                  } from '../modules/local/antismash'
include { DBCAN                                      } from '../modules/local/dbcan'
include { CIRCOS_PLOT                                } from '../modules/local/circos_plot'
include { PSEUDOFINDER                               } from '../modules/local/pseudofinder'
include { PSEUDOFINDER_POSTPROCESSING                } from '../modules/local/pseudofinder'

include { DOWNLOAD_DATABASES                         } from '../subworkflows/download_databases'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BAKTA_BAKTA                 } from '../modules/nf-core/bakta/bakta/main'
include { GECCO_RUN                   } from '../modules/nf-core/gecco/run/main'
include { QUAST                       } from '../modules/nf-core/quast/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

/////////////////////////////////////////////////////
/* --  Create channels for reference databases  -- */
/////////////////////////////////////////////////////

amrfinder_plus_db = channel.empty()

defense_finder_db = channel.empty()
dbcan_db = channel.empty()

interproscan_db = channel.empty()
interpro_entry_list = channel.empty()

eggnog_db = channel.empty()
eggnog_diamond_db = channel.empty()
eggnog_data = channel.empty()

rfam_ncrna_models = channel.empty()

pseudofinder_db = channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METTANNOTATOR {

    if (params.dbs) {
        // Download databases (if needed) //
        DOWNLOAD_DATABASES()

        amrfinder_plus_db = DOWNLOAD_DATABASES.out.amrfinder_plus_db

        antismash_db = DOWNLOAD_DATABASES.out.antismash_db

        defense_finder_db = DOWNLOAD_DATABASES.out.defense_finder_db

        dbcan_db = DOWNLOAD_DATABASES.out.dbcan_db

        interproscan_db = DOWNLOAD_DATABASES.out.interproscan_db

        interpro_entry_list = DOWNLOAD_DATABASES.out.interpro_entry_list

        eggnog_db = DOWNLOAD_DATABASES.out.eggnog_db

        rfam_ncrna_models = DOWNLOAD_DATABASES.out.rfam_ncrna_models

        pseudofinder_db = DOWNLOAD_DATABASES.out.pseudofinder_db

        if (params.bakta) {
            bakta_db = DOWNLOAD_DATABASES.out.bakta_db
        }
    } else {
        // Use the parametrized folders and files for the databases //
        amrfinder_plus_db = tuple(
            file(params.amrfinder_plus_db, checkIfExists: true),
            params.amrfinder_plus_db_version
        )

        antismash_db = tuple(
            file(params.antismash_db, checkIfExists: true),
            params.antismash_db_version
        )

        defense_finder_db = tuple(
            file(params.defense_finder_db, checkIfExists: true),
            params.defense_finder_db_version
        )

        dbcan_db = tuple(
            file(params.dbcan_db, checkIfExists: true),
            params.dbcan_db_version
        )

        interproscan_db = tuple(
            file(params.interproscan_db, checkIfExists: true),
            params.interproscan_db_version
        )

        interpro_entry_list = tuple(
            file(params.interpro_entry_list, checkIfExists: true),
            params.interpro_entry_list_version
        )

        eggnog_db = tuple(
            file(params.eggnog_db, checkIfExists: true),
            params.eggnog_db_version
        )

        rfam_ncrna_models = tuple(
            file(params.rfam_ncrna_models, checkIfExists: true),
            params.rfam_ncrna_models_rfam_version
        )

        pseudofinder_db = tuple(
            file(params.pseudofinder_db, checkIfExists: true),
            params.pseudofinder_db_version
        )

        if (params.bakta) {
            bakta_db = tuple(
                file(params.bakta_db, checkIfExists: true),
                params.bakta_db_version
            )
        }
    }

    ch_versions = Channel.empty()

    assemblies = Channel.fromSamplesheet("input")
    annotations_fna = channel.empty()
    annotations_gbk = channel.empty()
    annotations_faa = channel.empty()
    annotations_gff = channel.empty()
    compliant_gbk = channel.empty()
    compliant_gff = channel.empty()

    LOOKUP_KINGDOM( assemblies )
    ch_versions = ch_versions.mix(LOOKUP_KINGDOM.out.versions.first())

    assemblies_with_kingdom = assemblies.join( LOOKUP_KINGDOM.out.detected_kingdom )

   if ( params.bakta ) {
       assemblies_with_kingdom.branch {
           bacteria: it[2] == "Bacteria"
           archaea: it[2] == "Archaea"
       }.set { assemblies_to_annotate }

       BAKTA_BAKTA( assemblies_to_annotate.bacteria, bakta_db )

       PROKKA_STANDARD( assemblies_to_annotate.archaea, Channel.value("standard") )
       PROKKA_COMPLIANT( assemblies_to_annotate.archaea, Channel.value("compliant") )

       ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions.first())
       ch_versions = ch_versions.mix(PROKKA_STANDARD.out.versions.first())

       annotations_fna = annotations_fna.mix( BAKTA_BAKTA.out.fna ).mix( PROKKA_STANDARD.out.fna )
       annotations_gbk = annotations_gbk.mix( BAKTA_BAKTA.out.gbk ).mix( PROKKA_STANDARD.out.gbk )
       annotations_faa = annotations_faa.mix( BAKTA_BAKTA.out.faa ).mix( PROKKA_STANDARD.out.faa )
       annotations_gff = annotations_gff.mix( BAKTA_BAKTA.out.gff ).mix( PROKKA_STANDARD.out.gff )

       compliant_gbk = compliant_gbk.mix( BAKTA_BAKTA.out.gbk ).mix( PROKKA_COMPLIANT.out.gbk )
       compliant_gff = compliant_gff.mix( BAKTA_BAKTA.out.gff ).mix( PROKKA_COMPLIANT.out.gff )

   } else {

       PROKKA_STANDARD( assemblies_with_kingdom, Channel.value("standard") )
       PROKKA_COMPLIANT( assemblies_with_kingdom, Channel.value("compliant") )

       ch_versions = ch_versions.mix(PROKKA_STANDARD.out.versions.first())

       annotations_fna = PROKKA_STANDARD.out.fna
       annotations_gbk = PROKKA_STANDARD.out.gbk
       annotations_faa = PROKKA_STANDARD.out.faa
       annotations_gff = PROKKA_STANDARD.out.gff

       compliant_gbk = PROKKA_COMPLIANT.out.gbk
       compliant_gff = PROKKA_COMPLIANT.out.gff
   }

    assemblies_for_quast = assemblies.join(
        annotations_gff
    ).map { it -> tuple(it[0], it[1], it[2]) }

    QUAST(
        assemblies_for_quast.map { tuple(it[0], it[1]) }, // the assembly
        assemblies_for_quast.map { tuple(it[0], it[2]) }, // the GFF
    )

    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    PSEUDOFINDER(
        compliant_gbk,
        pseudofinder_db
    )

    ch_versions = ch_versions.mix(PSEUDOFINDER.out.versions.first())

    PSEUDOFINDER_POSTPROCESSING(
        annotations_gff.join(
            compliant_gff
        ).join(
            PSEUDOFINDER.out.pseudofinder_gff
        )
    )

    ch_versions = ch_versions.mix(PSEUDOFINDER_POSTPROCESSING.out.versions.first())

    CRISPRCAS_FINDER( assemblies )

    ch_versions = ch_versions.mix(CRISPRCAS_FINDER.out.versions.first())

    // EGGNOG_MAPPER_ORTHOLOGS - needs a third empty file in mode=mapper
    proteins_for_emapper_orth = annotations_faa.map { it -> tuple( it[0], file(it[1]), file("NO_FILE") ) }

    EGGNOG_MAPPER_ORTHOLOGS(
        proteins_for_emapper_orth,
        Channel.value("mapper"),
        eggnog_db
    )

    ch_versions = ch_versions.mix(EGGNOG_MAPPER_ORTHOLOGS.out.versions.first())

    // EGGNOG_MAPPER_ANNOTATIONS - needs a second empty file in mode=annotations
    orthologs_for_annotations = assemblies.join(EGGNOG_MAPPER_ORTHOLOGS.out.orthologs, by: 0).map {
        it -> {
            tuple(it[0], file("NO_FILE"), file(it[2])) // tuple( meta , <empty> , assembly )
        }
    }

    EGGNOG_MAPPER_ANNOTATIONS(
        orthologs_for_annotations,
        Channel.value("annotations"),
        eggnog_db
    )

    ch_versions = ch_versions.mix(EGGNOG_MAPPER_ANNOTATIONS.out.versions.first())

    if ( !params.fast ) {
        INTERPROSCAN(
            annotations_faa,
            interproscan_db
        )
        ch_versions = ch_versions.mix(INTERPROSCAN.out.versions.first())
    }

    assemblies_plus_faa_and_gff = assemblies.join(
        annotations_faa
    ).join(
        annotations_gff
    )

    AMRFINDER_PLUS(
        assemblies_plus_faa_and_gff,
        amrfinder_plus_db
    )

    ch_versions = ch_versions.mix(AMRFINDER_PLUS.out.versions.first())

    AMRFINDER_PLUS_TO_GFF( AMRFINDER_PLUS.out.amrfinder_tsv )

    ch_versions = ch_versions.mix(AMRFINDER_PLUS_TO_GFF.out.versions.first())

    DEFENSE_FINDER(
        annotations_faa.join( annotations_gff ),
        defense_finder_db
    )

    ch_versions = ch_versions.mix(DEFENSE_FINDER.out.versions.first())

    if ( !params.fast ) {
        UNIFIRE ( annotations_faa )
        ch_versions = ch_versions.mix(UNIFIRE.out.versions.first())
    }

    DETECT_TRNA(
        annotations_fna.join( LOOKUP_KINGDOM.out.detected_kingdom )
    )

    ch_versions = ch_versions.mix(DETECT_TRNA.out.versions.first())

    DETECT_NCRNA(
        annotations_fna,
        rfam_ncrna_models
    )

    ch_versions = ch_versions.mix(DETECT_NCRNA.out.versions.first())

    if ( !params.fast ) {
        SANNTIS(
            INTERPROSCAN.out.ips_annotations.join(annotations_gbk)
        )
        ch_versions = ch_versions.mix(SANNTIS.out.versions.first())
    }


    GECCO_RUN(
        annotations_gbk.map { meta, gbk -> [meta, gbk, []] }, []
    )

    ch_versions = ch_versions.mix(GECCO_RUN.out.versions.first())

    ANTISMASH(
        annotations_gbk,
        antismash_db
    )

    ch_versions = ch_versions.mix(ANTISMASH.out.versions.first())

    DBCAN(
        annotations_faa.join( annotations_gff ),
        dbcan_db
    )

    ch_versions = ch_versions.mix(DBCAN.out.versions.first())

    /**********************************************/
    /* Combine the results into a single GFF file */
    /**********************************************/
    annotate_gff_input = annotations_gff.join(
        EGGNOG_MAPPER_ANNOTATIONS.out.annotations
    ).join(
        DETECT_NCRNA.out.ncrna_tblout
    ).join(
        DETECT_TRNA.out.trna_gff
    ).join(
        CRISPRCAS_FINDER.out.hq_gff, remainder: true
    ).join(
        AMRFINDER_PLUS.out.amrfinder_tsv, remainder: true
    ).join(
        ANTISMASH.out.gff, remainder: true
    ).join(
        GECCO_RUN.out.gff, remainder: true
    ).join(
        DBCAN.out.dbcan_gff, remainder: true
    ).join(
        DEFENSE_FINDER.out.gff, remainder: true
    ).join(
        PSEUDOFINDER_POSTPROCESSING.out.pseudofinder_processed_gff, remainder: true
    )

    if ( !params.fast ) {
        annotate_gff_input = annotate_gff_input.join(
            INTERPROSCAN.out.ips_annotations
        ).join(
            SANNTIS.out.sanntis_gff, remainder: true
        ).join(
            UNIFIRE.out.arba, remainder: true
        ).join(
            UNIFIRE.out.unirule, remainder: true
        ).join(
            UNIFIRE.out.pirsr, remainder: true
        )
    } else {
        annotate_gff_input = annotate_gff_input.map { it -> {
                // IPS, SanntiS, UniFire{arba,unirule,pirsr}
                // meta, <files> //
                it + [[], [], [], [], []]
            }
        }
    }

    ANNOTATE_GFF(
        annotate_gff_input,
        interpro_entry_list
    )

    ch_versions = ch_versions.mix(ANNOTATE_GFF.out.versions.first())

    CIRCOS_PLOT(
        ANNOTATE_GFF.out.annotated_gff
    )

    ch_versions = ch_versions.mix(CIRCOS_PLOT.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMettannotator.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMettannotator.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix( QUAST.out.results.collect { it[1] }.ifEmpty([]) )


    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
