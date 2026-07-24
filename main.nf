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

    if (params.add_slow_tools && params.fast) {
        error "'--add_slow_tools' and '--fast' cannot be used together. Use '--fast' for initial run, then '--add-slow-tools' to add slow annotations."
    }

    if (params.skip_fast_staging && !params.add_slow_tools) {
        error "'--skip_fast_staging' can only be used together with '--add_slow_tools'."
    }

    if (params.allow_missing_files && !params.add_slow_tools) {
        error "'--allow_missing_files' can only be used together with '--add_slow_tools'."
    }

    if (params.ignore_version_mismatch && !params.add_slow_tools) {
        error "'--ignore_version_mismatch' can only be used together with '--add_slow_tools'."
    }

    if (params.gene_calling_only && params.fast) {
        error "'--gene_calling_only' and '--fast' cannot be used together."
    }

    if (params.gene_calling_only && params.add_slow_tools) {
        error "'--gene_calling_only' and '--add_slow_tools' cannot be used together."
    }

    if (params.gene_calling_only && params.gene_calls) {
        error "'--gene_calling_only' and '--gene_calls' cannot be used together."
    }

    if (params.gene_calls && params.add_slow_tools) {
        error "'--gene_calls' and '--add_slow_tools' cannot be used together."
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
