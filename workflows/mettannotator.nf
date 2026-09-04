/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { ANNOTATE_GFF                               } from '../modules/local/annotate_gff'
include { CIRCOS_PLOT                                } from '../modules/local/circos_plot'
include { STAGE_FAST_OUTPUTS                         } from '../modules/local/stage_fast_outputs'

include { DOWNLOAD_DATABASES                         } from '../subworkflows/download_databases'
include { EXTENDED_ANNOTATION                        } from '../subworkflows/extended_annotation'
include { GENE_CALLING                               } from '../subworkflows/gene_calling'
include { FAST_ANNOTATION                            } from '../subworkflows/fast_annotation'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow METTANNOTATOR {


    def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
    def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
    def summary_params = paramsSummaryMap(workflow)

    // Print parameter summary log to screen
    log.info logo + paramsSummaryLog(workflow) + citation

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

    bakta_db = channel.empty()

    // Info required for completion email and summary
    def multiqc_report = []


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        CONFIG FILES
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


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

    if ( params.extend_to_full ) {

        //
        // In extend-to-full mode, we read existing results from --fast and run only the additional tools
        //
       Channel.fromSamplesheet("input")
                    .map { row -> [row[0].prefix, row[0].taxid] }
                    .set { samplesheet }

        def results_dir = file(params.extend_to_full)

        def load_files = { pattern, strip_suffix ->
            channel
                .fromPath(pattern, checkIfExists: false)
                .map { f ->
                    def name = f.name.replaceAll(strip_suffix + '$', '')
                    tuple(name, f)
                }
        }

        // Load annotation files from either Prokka or Bakta output directories.
        // Prokka:  functional_annotation/prokka/  extensions .faa / .gbk / .gff
        // Bakta:   functional_annotation/bakta/   extensions .faa / .gbff / .gff3
        // Both channels are mixed so that a run using --bakta (or a mixed Prokka+Bakta
        // run, e.g. Bacteria→Bakta / Archaea→Prokka) is handled transparently.
        annotations_faa  = load_files("${results_dir}/**/functional_annotation/prokka/*.faa",  '\\.faa')
             .mix(load_files("${results_dir}/**/functional_annotation/bakta/*.faa",            '\\.faa'))
        annotations_gbk  = load_files("${results_dir}/**/functional_annotation/prokka/*.gbk",  '\\.gbk')
             .mix(load_files("${results_dir}/**/functional_annotation/bakta/*.gbff",           '\\.gbff'))
        annotations_gff  = load_files("${results_dir}/**/functional_annotation/prokka/*.gff",  '\\.gff')
             .mix(load_files("${results_dir}/**/functional_annotation/bakta/*.gff3",           '\\.gff3'))
        eggnog_files     = load_files("${results_dir}/**/functional_annotation/eggnog_mapper/*.emapper.annotations",'\\.emapper\\.annotations')
        ncrna_files      = load_files("${results_dir}/**/rnas/ncrna/*.ncrna.deoverlap.tbl",'\\.ncrna\\.deoverlap\\.tbl')
        trna_files       = load_files("${results_dir}/**/rnas/trna/*_trna.gff",'_trna\\.gff')
        crisprcas_files  = load_files("${results_dir}/**/mobilome/crisprcas_finder/*_crisprcasfinder_hq.gff",'_crisprcasfinder_hq\\.gff')
        amr_files        = load_files("${results_dir}/**/antimicrobial_resistance/amrfinder_plus/*_amrfinderplus.tsv",'_amrfinderplus\\.tsv')
        antismash_files  = load_files("${results_dir}/**/biosynthetic_gene_clusters/antismash/*_antismash.gff",'_antismash\\.gff')
        gecco_files      = load_files("${results_dir}/**/biosynthetic_gene_clusters/gecco/*_gecco_clusters.gff",'_gecco_clusters\\.gff')
        dbcan_files      = load_files("${results_dir}/**/functional_annotation/dbcan/*_dbcan.gff",'_dbcan\\.gff')
        df_files         = load_files("${results_dir}/**/antiphage_defense/defense_finder/*_defense_finder.gff",'_defense_finder\\.gff')
        pseudo_files     = load_files("${results_dir}/**/functional_annotation/pseudofinder/*_processed_pseudogenes.gff",'_processed_pseudogenes\\.gff')

        // Stage fast-run outputs into --outdir so the final directory
        // contains both fast and additional results in the expected layout.
        // Skipped when --skip-fast-staging is set.
        if ( !params.skip_fast_staging ) {
            stage_fast_input = samplesheet
                .map { prefix, taxid ->
                    def sample_dir = file("${results_dir}/${prefix}", checkIfExists: true)
                    tuple([prefix: prefix, taxid: taxid], sample_dir)
                }

            STAGE_FAST_OUTPUTS(stage_fast_input)
        }

        // Join samplesheet metadata with file paths
        samplesheet
            .join(annotations_faa)
            .map { prefix, taxid, faa ->
                def meta = [prefix: prefix, taxid: taxid]
                tuple(meta, faa)
            }
            .set { annotations_faa_input }
        samplesheet
            .join(annotations_gbk)
            .map { prefix, taxid, gbk ->
                def meta = [prefix: prefix, taxid: taxid]
                tuple(meta, gbk)
            }
            .set { annotations_gbk }

        EXTENDED_ANNOTATION(annotations_faa_input, annotations_gbk, interproscan_db)
        ch_versions = ch_versions.mix(EXTENDED_ANNOTATION.out.versions)

        // Re-key additional-tool outputs by prefix string for joining with disk-loaded fast results
        ips_by_prefix     = EXTENDED_ANNOTATION.out.ips_annotations.map { meta, f -> [meta.prefix, f] }
        sanntis_by_prefix = EXTENDED_ANNOTATION.out.sanntis_gff.map     { meta, f -> [meta.prefix, f] }
        arba_by_prefix    = EXTENDED_ANNOTATION.out.arba.map            { meta, f -> [meta.prefix, f] }
        unirule_by_prefix = EXTENDED_ANNOTATION.out.unirule.map         { meta, f -> [meta.prefix, f] }
        pirsr_by_prefix   = EXTENDED_ANNOTATION.out.pirsr.map           { meta, f -> [meta.prefix, f] }

        // Build fast results channel keyed by prefix
        fast_results = samplesheet
            .join(annotations_gff)
            .join(eggnog_files,    remainder: true)
            .join(ncrna_files,     remainder: true)
            .join(trna_files,      remainder: true)
            .join(crisprcas_files, remainder: true)
            .join(amr_files,       remainder: true)
            .join(antismash_files, remainder: true)
            .join(gecco_files,     remainder: true)
            .join(dbcan_files,     remainder: true)
            .join(df_files,        remainder: true)
            .join(pseudo_files,    remainder: true)

        // Merge fast + extended, then build final meta tuple
        annotate_gff_input = fast_results
            .join(ips_by_prefix)
            .join(sanntis_by_prefix, remainder: true)
            .join(arba_by_prefix,    remainder: true)
            .join(unirule_by_prefix, remainder: true)
            .join(pirsr_by_prefix,   remainder: true)
            .map { prefix, taxid,
                gff, eggnog, ncrna, trna, crisprcas, amr,
                antismash, gecco, dbcan, df, pseudo,
                ips, sanntis, arba, unirule, pirsr ->
                tuple(
                    [prefix: prefix, taxid: taxid],
                    gff, eggnog, ncrna, trna, crisprcas, amr,
                    antismash, gecco, dbcan, df, pseudo,
                    ips, sanntis, arba, unirule, pirsr
                )
            }

        ANNOTATE_GFF(
            annotate_gff_input,
            interpro_entry_list
        )
        ch_versions = ch_versions.mix(ANNOTATE_GFF.out.versions.first())

        CUSTOM_DUMPSOFTWAREVERSIONS(
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

    } else {

        assemblies_raw = Channel.fromSamplesheet("input")
        assemblies = assemblies_raw.map { row -> [row[0], row[1]] }

        GENE_CALLING(assemblies_raw, assemblies, bakta_db)
        ch_versions = ch_versions.mix(GENE_CALLING.out.versions)

        annotations_faa  = GENE_CALLING.out.annotations_faa
        annotations_gff  = GENE_CALLING.out.annotations_gff
        annotations_gbk  = GENE_CALLING.out.annotations_gbk
        annotations_fna  = GENE_CALLING.out.annotations_fna
        compliant_gbk    = GENE_CALLING.out.compliant_gbk
        compliant_gff    = GENE_CALLING.out.compliant_gff

        if ( !params.gene_calling_only ) {

            FAST_ANNOTATION(
                annotations_faa,
                annotations_gff,
                annotations_gbk,
                annotations_fna,
                compliant_gbk,
                compliant_gff,
                assemblies,
                GENE_CALLING.out.detected_kingdom,
                amrfinder_plus_db,
                antismash_db,
                defense_finder_db,
                dbcan_db,
                eggnog_db,
                rfam_ncrna_models,
                pseudofinder_db
            )
            ch_versions = ch_versions.mix(FAST_ANNOTATION.out.versions)

            if ( !params.fast ) {
                EXTENDED_ANNOTATION(annotations_faa, annotations_gbk, interproscan_db)
                ch_versions = ch_versions.mix(EXTENDED_ANNOTATION.out.versions)
            }

            /**********************************************/
            /* Combine the results into a single GFF file */
            /**********************************************/
            annotate_gff_input = annotations_gff.join(
                FAST_ANNOTATION.out.eggnog_annotations
            ).join(
                FAST_ANNOTATION.out.ncrna_tblout
            ).join(
                FAST_ANNOTATION.out.trna_gff
            ).join(
                FAST_ANNOTATION.out.crisprcas_hq_gff, remainder: true
            ).join(
                FAST_ANNOTATION.out.amrfinder_tsv, remainder: true
            ).join(
                FAST_ANNOTATION.out.antismash_gff, remainder: true
            ).join(
                FAST_ANNOTATION.out.gecco_gff, remainder: true
            ).join(
                FAST_ANNOTATION.out.dbcan_gff, remainder: true
            ).join(
                FAST_ANNOTATION.out.defense_gff, remainder: true
            ).join(
                FAST_ANNOTATION.out.pseudofinder_gff, remainder: true
            )

            if ( !params.fast ) {
                annotate_gff_input = annotate_gff_input.join(
                    EXTENDED_ANNOTATION.out.ips_annotations
                ).join(
                    EXTENDED_ANNOTATION.out.sanntis_gff, remainder: true
                ).join(
                    EXTENDED_ANNOTATION.out.arba, remainder: true
                ).join(
                    EXTENDED_ANNOTATION.out.unirule, remainder: true
                ).join(
                    EXTENDED_ANNOTATION.out.pirsr, remainder: true
                )
            } else {
                annotate_gff_input = annotate_gff_input.map { it ->
                        // IPS, SanntiS, UniFire{arba,unirule,pirsr}
                        // meta, <files> //
                        it + [[], [], [], [], []]
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

        }

        CUSTOM_DUMPSOFTWAREVERSIONS(
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

        if ( !params.gene_calling_only ) {
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
            ch_multiqc_files = ch_multiqc_files.mix( FAST_ANNOTATION.out.quast_results.collect { it[1] }.ifEmpty([]) )


            MULTIQC(
                ch_multiqc_files.collect(),
                ch_multiqc_config.toList(),
                ch_multiqc_custom_config.toList(),
                ch_multiqc_logo.toList()
            )
            multiqc_report = MULTIQC.out.report.toList()
        }
    }
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
