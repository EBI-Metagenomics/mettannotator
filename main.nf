#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ebi-metagenomics/mettannotator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ebi-metagenomics/mettannotator
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METTANNOTATOR } from './workflows/mettannotator'

//
// WORKFLOW: Run main ebi-metagenomics/mettannotator analysis pipeline
//
workflow EBIMETAGENOMICS_METTANNOTATOR {

    // Print help message if needed
    if (params.help) {
        def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
        def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
        def String command = "nextflow run ${workflow.manifest.name} --input assemblies_sheet.csv -profile docker"
        log.info logo + paramsHelp(command) + citation + NfcoreTemplate.dashedLine(params.monochrome_logs)
        System.exit(0)
    }

    // Validate input parameters
    if (params.validate_params) {
        validateParameters()
    }

    // Custom validation until conditional schema validation gets implemented
    // https://github.com/nf-core/tools/issues/2453
    if (params.dbs == null && (
        params.amrfinder_plus_db == null ||
        params.antismash_db == null ||
        params.defense_finder_db == null ||
        params.dbcan_db == null ||
        params.interproscan_db == null ||
        params.interpro_entry_list == null ||
        params.eggnog_db == null ||
        params.rfam_ncrna_models == null
    )) {
        error "If the parameter '--dbs' is null, you must specify individual paths for each database."
    }

    if (params.extend_to_full && params.fast) {
        error "'--extend_to_full' and '--fast' cannot be used together. Use '--fast' for initial run, then '--extend_to_full' to add extended annotations."
    }

    if (params.skip_fast_staging && !params.extend_to_full) {
        error "'--skip_fast_staging' can only be used together with '--extend_to_full'."
    }

    if (params.allow_missing_files && !params.extend_to_full) {
        error "'--allow_missing_files' can only be used together with '--extend_to_full'."
    }

    if (params.ignore_version_mismatch && !params.extend_to_full) {
        error "'--ignore_version_mismatch' can only be used together with '--extend_to_full'."
    }

    if (params.gene_calling_only && params.fast) {
        error "'--gene_calling_only' and '--fast' cannot be used together."
    }

    if (params.gene_calling_only && params.extend_to_full) {
        error "'--gene_calling_only' and '--extend_to_full' cannot be used together."
    }

    if (params.gene_calling_only && params.gene_calls) {
        error "'--gene_calling_only' and '--gene_calls' cannot be used together."
    }

    if (params.gene_calls && params.extend_to_full) {
        error "'--gene_calls' and '--extend_to_full' cannot be used together."
    }

    if (params.gene_calls && params.input) {
        def ssHeader = new File(params.input as String).readLines()[0].split(',')*.trim()
        if (ssHeader.contains('gene_calls_gff')) {
           error "'--gene_calls' cannot be used when the samplesheet contains a 'gene_calls_gff' column. Use one mechanism or the other, not both."
        }
    }

    WorkflowMain.initialise(workflow, params, log)
    WorkflowMettannotator.initialise(params, log, workflow)

    METTANNOTATOR ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    EBIMETAGENOMICS_METTANNOTATOR ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
