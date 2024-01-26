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
include { LOOKUP_KINGDOM } from '../modules/local/lookup_kingdom'
include { PROKKA } from '../modules/local/prokka'
include { AMRFINDER_PLUS; AMRFINDER_PLUS_TO_GFF } from '../modules/local/amrfinder_plus'
include { DEFENSE_FINDER } from '../modules/local/defense_finder'
include { CRISPRCAS_FINDER } from '../modules/local/crisprcasfinder'
include { EGGNOG_MAPPER as EGGNOG_MAPPER_ORTHOLOGS } from '../modules/local/eggnog'
include { EGGNOG_MAPPER as EGGNOG_MAPPER_ANNOTATIONS } from '../modules/local/eggnog'
include { INTERPROSCAN } from '../modules/local/interproscan'
include { DETECT_RRNA } from '../modules/local/detect_rrna'
include { DETECT_NCRNA } from '../modules/local/detect_ncrna'
include { SANNTIS } from '../modules/local/sanntis'
include { UNIFIRE } from '../modules/local/unifire'
include { ANNOTATE_GFF } from '../modules/local/annotate_gff'
include { ANTISMASH } from '../modules/local/antismash'
include { DBCAN } from '../modules/local/dbcan'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GECCO_RUN } from '../modules/nf-core/gecco/run/main'
include { QUAST } from '../modules/nf-core/quast/main'
include { MULTIQC } from '../modules/nf-core/multiqc/main'
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

ch_interproscan_db = file(params.interproscan_db)

ch_eggnog_db = file(params.eggnog_db)
ch_eggnog_diamond_db = file(params.eggnong_diamond_db)
ch_eggnog_data_dir = file(params.eggnong_data_dir)

ch_rfam_rrna_models = file(params.rfam_rrna_models)
ch_rfam_ncrna_models = file(params.rfam_ncrna_models)

ch_amrfinder_plus_db = file(params.amrfinder_plus_db)

ch_defense_finder_db = file(params.defense_finder_db)

ch_dbcan_db = file(params.dbcan_db)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METTANNOTATOR {

    ch_versions = Channel.empty()

    assemblies = Channel.fromSamplesheet("input")

    LOOKUP_KINGDOM ( )

    PROKKA( assemblies.join( LOOKUP_KINGDOM.detected_kingdom ))

    ch_versions = ch_versions.mix(PROKKA.out.versions.first())

    assemblies_for_quast = assemblies.join(
        PROKKA.out.gff
    ).map { it -> tuple(it[0], it[1], it[2]) }

    QUAST(
        assemblies_for_quast.map { tuple(it[0], it[1]) }, // the assembly
        assemblies_for_quast.map { tuple(it[0], it[2]) }, // the GFF
    )

    ch_versions = ch_versions.mix(QUAST.out.versions.first())

    CRISPRCAS_FINDER( assemblies )

    ch_versions = ch_versions.mix(CRISPRCAS_FINDER.out.versions.first())

    // EGGNOG_MAPPER_ORTHOLOGS - needs a third empty file in mode=mapper
    proteins_for_emapper_orth = PROKKA.out.faa.map { it -> tuple( it[0], file(it[1]), file("NO_FILE") ) }

    EGGNOG_MAPPER_ORTHOLOGS(
        proteins_for_emapper_orth,
        Channel.value("mapper"),
        ch_eggnog_db,
        ch_eggnog_diamond_db,
        ch_eggnog_data_dir
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
        ch_eggnog_db,
        ch_eggnog_diamond_db,
        ch_eggnog_data_dir
    )

    ch_versions = ch_versions.mix(EGGNOG_MAPPER_ANNOTATIONS.out.versions.first())

    INTERPROSCAN(
        PROKKA.out.faa,
        ch_interproscan_db
    )

    ch_versions = ch_versions.mix(INTERPROSCAN.out.versions.first())

    assemblies_plus_faa_and_gff = assemblies.join(
        PROKKA.out.faa
    ).join(
        PROKKA.out.gff
    )

    AMRFINDER_PLUS( assemblies_plus_faa_and_gff )

    ch_versions = ch_versions.mix(AMRFINDER_PLUS.out.versions.first())

    AMRFINDER_PLUS_TO_GFF( AMRFINDER_PLUS.out.amrfinder_tsv )

    ch_versions = ch_versions.mix(AMRFINDER_PLUS_TO_GFF.out.versions.first())

    DEFENSE_FINDER (
        PROKKA.out.faa.join( PROKKA.out.gff ),
        ch_defense_finder_db
    )

    ch_versions = ch_versions.mix(DEFENSE_FINDER.out.versions.first())

    UNIFIRE ( PROKKA.out.faa )

    ch_versions = ch_versions.mix(UNIFIRE.out.versions.first())

    DETECT_RRNA(
        PROKKA.out.fna,
        ch_rfam_rrna_models
    )

    ch_versions = ch_versions.mix(DETECT_RRNA.out.versions.first())

    DETECT_NCRNA(
        PROKKA.out.fna,
        ch_rfam_ncrna_models
    )

    ch_versions = ch_versions.mix(DETECT_NCRNA.out.versions.first())

    SANNTIS(
        INTERPROSCAN.out.ips_annotations.join(PROKKA.out.gbk)
    )

    ch_versions = ch_versions.mix(SANNTIS.out.versions.first())

    GECCO_RUN(
        PROKKA.out.gbk.map { meta, gbk -> [meta, gbk, []] }, []
    )

    ch_versions = ch_versions.mix(GECCO_RUN.out.versions.first())

    ANTISMASH(
        PROKKA.out.gbk
    )

    ch_versions = ch_versions.mix(ANTISMASH.out.versions.first())

    DBCAN(
        PROKKA.out.faa.join( PROKKA.out.gff ),
        ch_dbcan_db
    )

    ch_versions = ch_versions.mix(DBCAN.out.versions.first())

    ANNOTATE_GFF(
        PROKKA.out.gff.join(
            INTERPROSCAN.out.ips_annotations
        ).join(
           EGGNOG_MAPPER_ANNOTATIONS.out.annotations, remainder: true
        ).join(
            SANNTIS.out.sanntis_gff, remainder: true
        ).join(
            DETECT_NCRNA.out.ncrna_tblout
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
            UNIFIRE.out.arba, remainder: true
        ).join(
            UNIFIRE.out.unirule, remainder: true
        ).join(
            UNIFIRE.out.pirsr, remainder: true
        )
    )

    ch_versions = ch_versions.mix(ANNOTATE_GFF.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
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
    ch_multiqc_files = ch_multiqc_files.mix( PROKKA.out.txt.collect { it[1] }.ifEmpty([]) )

    MULTIQC (
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
