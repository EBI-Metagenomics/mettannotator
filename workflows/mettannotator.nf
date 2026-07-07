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
include { LOOKUP_KINGDOM                             } from '../modules/local/lookup_kingdom'
include { PROKKA as PROKKA_STANDARD                  } from '../modules/local/prokka'
include { PROKKA as PROKKA_COMPLIANT                 } from '../modules/local/prokka'
include { RENAME_CONTIGS as SHORTEN_CONTIG_NAMES     } from '../modules/local/rename_contigs'
include { RENAME_CONTIGS as REVERT_CONTIG_RENAMING   } from '../modules/local/rename_contigs'
include { AMRFINDER_PLUS_GET_SPECIES_NAME            } from '../modules/local/amrfinder_plus'
include { AMRFINDER_PLUS; AMRFINDER_PLUS_TO_GFF      } from '../modules/local/amrfinder_plus'
include { AMRFINDER_PLUS_TSV_POSTPROCESSING          } from '../modules/local/amrfinder_plus'
include { DEFENSE_FINDER                             } from '../modules/local/defense_finder'
include { CRISPRCAS_FINDER                           } from '../modules/local/crisprcasfinder'
include { EGGNOG_MAPPER as EGGNOG_MAPPER_ORTHOLOGS   } from '../modules/local/eggnog'
include { EGGNOG_MAPPER as EGGNOG_MAPPER_ANNOTATIONS } from '../modules/local/eggnog'
include { ADD_TAXID_TO_PROTEIN_FASTA                 } from '../modules/local/add_taxid'
include { INTERPROSCAN                               } from '../modules/local/interproscan'
include { DETECT_TRNA                                } from '../modules/local/detect_trna'
include { DETECT_NCRNA                               } from '../modules/local/detect_ncrna'
include { SANNTIS                                    } from '../modules/local/sanntis'
include { UNIFIRE                                    } from '../modules/local/unifire'
include { ANNOTATE_GFF                               } from '../modules/local/annotate_gff'
include { ANTISMASH                                  } from '../modules/local/antismash'
include { ANTISMASH_TO_GFF                           } from '../modules/local/antismash'
include { ANTISMASH_SUMMARY                          } from '../modules/local/antismash'
include { DBCAN                                      } from '../modules/local/dbcan'
include { CIRCOS_PLOT                                } from '../modules/local/circos_plot'
include { PSEUDOFINDER                               } from '../modules/local/pseudofinder'
include { PSEUDOFINDER_POSTPROCESSING                } from '../modules/local/pseudofinder'
include { STAGE_FAST_OUTPUTS                         } from '../modules/local/stage_fast_outputs'
include { NORMALIZE_ENSEMBL_GFF                      } from '../modules/local/normalize_ensembl_gff'
include { GFF_TO_GBK                                 } from '../modules/local/gff_to_gbk'


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

    if ( params.add_slow_tools ) {

        //
        // In add-slow-tools mode, we read existing results from --fast and run only slow tools
        //
        Channel.fromSamplesheet("input")
                     .map { meta, assembly ->
                        [meta.prefix, meta.taxid]
                    }
                    .set { samplesheet }

        def results_dir = file(params.add_slow_tools)

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
        // contains both fast and slow results in the expected layout.
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

        // For InterProScan, we need to add taxid to the protein FASTA
        ADD_TAXID_TO_PROTEIN_FASTA(
            annotations_faa_input
        )
        ch_versions = ch_versions.mix(ADD_TAXID_TO_PROTEIN_FASTA.out.versions.first())

        INTERPROSCAN(
            ADD_TAXID_TO_PROTEIN_FASTA.out.annotations_faa_with_taxid,
            interproscan_db
        )
        ch_versions = ch_versions.mix(INTERPROSCAN.out.versions.first())

        UNIFIRE (
            INTERPROSCAN.out.ips_xml
        )
        ch_versions = ch_versions.mix(UNIFIRE.out.versions.first())

        SANNTIS(
           INTERPROSCAN.out.ips_annotations.join(annotations_gbk)
        )
        ch_versions = ch_versions.mix(SANNTIS.out.versions.first())

        // Re-key slow tool outputs by prefix string for joining
        ips_by_prefix     = INTERPROSCAN.out.ips_annotations.map { meta, f -> [meta.prefix, f] }
        sanntis_by_prefix = SANNTIS.out.sanntis_gff.map         { meta, f -> [meta.prefix, f] }
        arba_by_prefix    = UNIFIRE.out.arba.map                { meta, f -> [meta.prefix, f] }
        unirule_by_prefix = UNIFIRE.out.unirule.map             { meta, f -> [meta.prefix, f] }
        pirsr_by_prefix   = UNIFIRE.out.pirsr.map               { meta, f -> [meta.prefix, f] }

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

        // Merge fast + slow, then build final meta tuple
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
        // [meta, assembly] for all samples — feeds LOOKUP_KINGDOM and all assembly-level tools
        assemblies = assemblies_raw.map { row -> [row[0], row[1]] }
        annotations_fna = channel.empty()
        annotations_gbk = channel.empty()
        annotations_faa = channel.empty()
        annotations_gff = channel.empty()
        compliant_gbk = channel.empty()
        compliant_gff = channel.empty()

        // Per-sample external gene calls: rows where annotation_source is prokka/bakta/ensembl.
        // When the gene_calls_gff column is absent (standard samplesheet), this channel is empty.
        // Tuple: [meta, assembly, gene_calls_gff, gene_calls_faa, gene_calls_gbk]
        ext_gene_call_rows = assemblies_raw
            .filter { it[0].annotation_source in ['prokka', 'bakta', 'ensembl'] }
            .map { [it[0], it[1], it[2], it[3], it[4]] }

        ext_gene_call_rows.branch {
            ensembl: it[0].annotation_source == 'ensembl'
            direct:  true   // prokka or bakta — GFF used as-is
        }.set { ext_gc_branches }

        NORMALIZE_ENSEMBL_GFF(
            ext_gc_branches.ensembl.map { meta, _assembly, gff, _faa, _gbk -> [meta, gff] }
        )
        ch_versions = ch_versions.mix(NORMALIZE_ENSEMBL_GFF.out.versions)

        // Ensembl samples with a user-supplied GBK skip GFF_TO_GBK.
        ext_gc_branches.ensembl.branch {
            has_gbk:   it[4]
            needs_gbk: true
        }.set { ensembl_gbk_branches }

        GFF_TO_GBK(
            NORMALIZE_ENSEMBL_GFF.out.normalised_gff
                .join( ensembl_gbk_branches.needs_gbk.map { meta, _assembly, _gff, faa, _gbk -> [meta, faa] } )
                .join( ensembl_gbk_branches.needs_gbk.map { meta, assembly, _gff, _faa, _gbk -> [meta, assembly] } )
        )
        ch_versions = ch_versions.mix(GFF_TO_GBK.out.versions)

        ext_annotations_gff = NORMALIZE_ENSEMBL_GFF.out.normalised_gff
            .mix( ext_gc_branches.direct.map { meta, _assembly, gff, _faa, _gbk -> [meta, gff] } )
        ext_annotations_faa = ext_gc_branches.ensembl.map { meta, _assembly, _gff, faa, _gbk -> [meta, faa] }
            .mix( ext_gc_branches.direct.map { meta, _assembly, _gff, faa, _gbk -> [meta, faa] } )
        ext_annotations_fna = ext_gc_branches.ensembl.map { meta, assembly, _gff, _faa, _gbk -> [meta, assembly] }
            .mix( ext_gc_branches.direct.map { meta, assembly, _gff, _faa, _gbk -> [meta, assembly] } )

        // Any external sample (ensembl or direct) that provides a GBK uses it directly,
        // bypassing GFF_TO_GBK. compliant_gff for direct samples uses their GFF as-is
        // (no contig renaming). Ensembl compliant_gff comes from NORMALIZE_ENSEMBL_GFF.
        ext_gbk = ensembl_gbk_branches.has_gbk
            .map    { meta, _assembly, _gff, _faa, gbk -> [meta, gbk] }
            .mix( ext_gc_branches.direct
                .filter { meta, _assembly, _gff, _faa, gbk -> gbk }
                .map    { meta, _assembly, _gff, _faa, gbk -> [meta, gbk] } )
        ext_direct_compliant_gff = ext_gc_branches.direct
            .filter { meta, _assembly, _gff, _faa, gbk -> gbk }
            .map    { meta, _assembly, gff, _faa, _gbk -> [meta, gff] }

        LOOKUP_KINGDOM( assemblies )
        ch_versions = ch_versions.mix(LOOKUP_KINGDOM.out.versions.first())

        assemblies_with_kingdom = assemblies.join( LOOKUP_KINGDOM.out.detected_kingdom )

        if ( params.gene_calls ) {

            def gene_calls_dir = file(params.gene_calls)

            def load_files = { pattern, strip_suffix ->
                channel
                    .fromPath(pattern, checkIfExists: false)
                    .map { f ->
                        def name = f.name.replaceAll(strip_suffix + '$', '')
                        tuple(name, f)
                    }
            }

            // Preserve the original samplesheet meta (needed for downstream joins, e.g. LOOKUP_KINGDOM)
            assemblies_by_prefix = assemblies.map { meta, assembly -> [meta.prefix, meta] }

            // Load annotation files from the appropriate tool's output directory.
            // When --bakta is set, load from functional_annotation/bakta/ (bakta compliant = standard bakta outputs).
            // Otherwise load from functional_annotation/prokka/ and functional_annotation/prokka_compliant/.
            // Prokka: .faa / .gbk / .gff / .fna   Bakta: .faa / .gbff / .gff3 / .fna
            if ( params.bakta ) {
                gc_faa           = load_files("${gene_calls_dir}/**/functional_annotation/bakta/*.faa",   '\\.faa')
                gc_gbk           = load_files("${gene_calls_dir}/**/functional_annotation/bakta/*.gbff",  '\\.gbff')
                gc_gff           = load_files("${gene_calls_dir}/**/functional_annotation/bakta/*.gff3",  '\\.gff3')
                gc_fna           = load_files("${gene_calls_dir}/**/functional_annotation/bakta/*.fna",   '\\.fna')
                gc_compliant_gbk = load_files("${gene_calls_dir}/**/functional_annotation/bakta/*.gbff",  '\\.gbff')
                gc_compliant_gff = load_files("${gene_calls_dir}/**/functional_annotation/bakta/*.gff3",  '\\.gff3')
            } else {
                gc_faa           = load_files("${gene_calls_dir}/**/functional_annotation/prokka/*.faa",            '\\.faa')
                gc_gbk           = load_files("${gene_calls_dir}/**/functional_annotation/prokka/*.gbk",            '\\.gbk')
                gc_gff           = load_files("${gene_calls_dir}/**/functional_annotation/prokka/*.gff",            '\\.gff')
                gc_fna           = load_files("${gene_calls_dir}/**/functional_annotation/prokka/*.fna",            '\\.fna')
                gc_compliant_gbk = load_files("${gene_calls_dir}/**/functional_annotation/prokka_compliant/*.gbk", '\\.gbk')
                gc_compliant_gff = load_files("${gene_calls_dir}/**/functional_annotation/prokka_compliant/*.gff", '\\.gff')
            }

            annotations_faa  = assemblies_by_prefix.join(gc_faa).map           { prefix, meta, f -> tuple(meta, f) }
            annotations_gbk  = assemblies_by_prefix.join(gc_gbk).map           { prefix, meta, f -> tuple(meta, f) }
            annotations_gff  = assemblies_by_prefix.join(gc_gff).map           { prefix, meta, f -> tuple(meta, f) }
            annotations_fna  = assemblies_by_prefix.join(gc_fna).map           { prefix, meta, f -> tuple(meta, f) }
            compliant_gbk    = assemblies_by_prefix.join(gc_compliant_gbk).map { prefix, meta, f -> tuple(meta, f) }
            compliant_gff    = assemblies_by_prefix.join(gc_compliant_gff).map { prefix, meta, f -> tuple(meta, f) }

        } else {

            // Exclude samples that supply per-sample external gene calls from Prokka/Bakta
            assemblies_for_gene_calling = assemblies_with_kingdom.filter { meta, _assembly, _kingdom ->
                !(meta.annotation_source in ['prokka', 'bakta', 'ensembl'])
            }

            if ( params.bakta ) {

                assemblies_for_gene_calling.branch {
                    bacteria: it[2] == "Bacteria"
                    archaea: it[2] == "Archaea"
                }.set { assemblies_to_annotate }

                BAKTA_BAKTA( assemblies_to_annotate.bacteria, bakta_db )
                ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions.first())

                prokka_input = assemblies_to_annotate.archaea

            } else {

                prokka_input = assemblies_for_gene_calling

            }

            // Prokka crashes with long contig names, so we temporary rename contigs before running it
            SHORTEN_CONTIG_NAMES(
                prokka_input.map { meta, fasta, _kingdom -> [ meta, fasta ] },
                []
            )
            ch_versions = ch_versions.mix(SHORTEN_CONTIG_NAMES.out.versions.first())

            renamed_prokka_input = SHORTEN_CONTIG_NAMES.out.modified_fasta
                .join( prokka_input )
                .map { meta, renamed_fasta, _original_fasta, kingdom ->
                    [ meta, renamed_fasta, kingdom ]
                }
            PROKKA_STANDARD( renamed_prokka_input, Channel.value("standard") )
            PROKKA_COMPLIANT( renamed_prokka_input, Channel.value("compliant") )
            ch_versions = ch_versions.mix(PROKKA_STANDARD.out.versions.first())

            // Revert contig renaming in both GBK and GFF files, so they match the original FASTA
            PROKKA_STANDARD.out.gbk
                .mix(PROKKA_STANDARD.out.gff)
                .mix(PROKKA_STANDARD.out.fna)
                .groupTuple()
                .join(SHORTEN_CONTIG_NAMES.out.fasta_ids_mapping)
                .multiMap { meta, files_list, names_mapping ->
                    files: [ meta, files_list ]
                    mapping: names_mapping
            }.set { revert_contig_renaming_input }
            REVERT_CONTIG_RENAMING(
                revert_contig_renaming_input.files,
                revert_contig_renaming_input.mapping
            )

            if ( params.bakta ) {

                annotations_gbk = annotations_gbk.mix( BAKTA_BAKTA.out.gbk ).mix( REVERT_CONTIG_RENAMING.out.modified_gbk )
                annotations_gff = annotations_gff.mix( BAKTA_BAKTA.out.gff ).mix( REVERT_CONTIG_RENAMING.out.modified_gff )
                annotations_fna = annotations_fna.mix( BAKTA_BAKTA.out.fna ).mix( REVERT_CONTIG_RENAMING.out.modified_fasta )
                annotations_faa = annotations_faa.mix( BAKTA_BAKTA.out.faa ).mix( PROKKA_STANDARD.out.faa )
                compliant_gbk = compliant_gbk.mix( BAKTA_BAKTA.out.gbk ).mix( PROKKA_COMPLIANT.out.gbk )
                compliant_gff = compliant_gff.mix( BAKTA_BAKTA.out.gff ).mix( PROKKA_COMPLIANT.out.gff )

            } else {

                annotations_gbk = REVERT_CONTIG_RENAMING.out.modified_gbk
                annotations_gff = REVERT_CONTIG_RENAMING.out.modified_gff
                annotations_fna = REVERT_CONTIG_RENAMING.out.modified_fasta
                annotations_faa = PROKKA_STANDARD.out.faa

                compliant_gbk = PROKKA_COMPLIANT.out.gbk
                compliant_gff = PROKKA_COMPLIANT.out.gff
            }
        }

        // Mix per-sample external gene calls (samplesheet gene_calls_gff column) into the main
        // annotation channels. Empty no-op when no external gene calls are provided.
        annotations_gff = annotations_gff.mix(ext_annotations_gff)
        annotations_faa = annotations_faa.mix(ext_annotations_faa)
        annotations_fna = annotations_fna.mix(ext_annotations_fna)
        annotations_gbk = annotations_gbk.mix(GFF_TO_GBK.out.gbk).mix(ext_gbk)
        compliant_gbk   = compliant_gbk.mix(GFF_TO_GBK.out.gbk).mix(ext_gbk)
        compliant_gff   = compliant_gff.mix(NORMALIZE_ENSEMBL_GFF.out.normalised_gff).mix(ext_direct_compliant_gff)


        if ( !params.gene_calling_only ) {

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
                ADD_TAXID_TO_PROTEIN_FASTA(
                    annotations_faa
                )
                ch_versions = ch_versions.mix(ADD_TAXID_TO_PROTEIN_FASTA.out.versions.first())

                INTERPROSCAN(
                    ADD_TAXID_TO_PROTEIN_FASTA.out.annotations_faa_with_taxid,
                    interproscan_db
                )
                ch_versions = ch_versions.mix(INTERPROSCAN.out.versions.first())

                UNIFIRE (
                    INTERPROSCAN.out.ips_xml
                )
                ch_versions = ch_versions.mix(UNIFIRE.out.versions.first())
            }

            assemblies_plus_faa_and_gff = assemblies.join(
                annotations_faa
            ).join(
                annotations_gff
            )

            AMRFINDER_PLUS_GET_SPECIES_NAME(
                assemblies,
                amrfinder_plus_db
            )

            AMRFINDER_PLUS(
                assemblies_plus_faa_and_gff.join(AMRFINDER_PLUS_GET_SPECIES_NAME.out.detected_organism),
                amrfinder_plus_db
            )
            ch_versions = ch_versions.mix(AMRFINDER_PLUS.out.versions.first())

            AMRFINDER_PLUS_TSV_POSTPROCESSING(AMRFINDER_PLUS.out.amrfinder_tsv)
            ch_versions = ch_versions.mix(AMRFINDER_PLUS_TSV_POSTPROCESSING.out.versions.first())

            AMRFINDER_PLUS_TO_GFF( AMRFINDER_PLUS_TSV_POSTPROCESSING.out.amrfinder_tsv )
            ch_versions = ch_versions.mix(AMRFINDER_PLUS_TO_GFF.out.versions.first())

            DEFENSE_FINDER(
                annotations_faa.join( annotations_gff ),
                defense_finder_db
            )
            ch_versions = ch_versions.mix(DEFENSE_FINDER.out.versions.first())

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

            ANTISMASH_TO_GFF(
                ANTISMASH.out.json
            )

            ch_versions = ch_versions.mix(ANTISMASH_TO_GFF.out.versions.first())

            ANTISMASH_SUMMARY(
                ANTISMASH_TO_GFF.out.gff
            )
            ch_versions = ch_versions.mix(ANTISMASH_SUMMARY.out.versions.first())

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
                ANTISMASH_TO_GFF.out.gff, remainder: true
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
            ch_multiqc_files = ch_multiqc_files.mix( QUAST.out.results.collect { it[1] }.ifEmpty([]) )


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
