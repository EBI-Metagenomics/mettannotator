include { LOOKUP_KINGDOM                             } from '../modules/local/lookup_kingdom'
include { PROKKA as PROKKA_STANDARD                  } from '../modules/local/prokka'
include { PROKKA as PROKKA_COMPLIANT                 } from '../modules/local/prokka'
include { RENAME_CONTIGS as SHORTEN_CONTIG_NAMES     } from '../modules/local/rename_contigs'
include { RENAME_CONTIGS as REVERT_CONTIG_RENAMING   } from '../modules/local/rename_contigs'
include { BAKTA_BAKTA                                } from '../modules/nf-core/bakta/bakta/main'
include { NORMALIZE_ENSEMBL_GFF                      } from '../modules/local/normalize_ensembl_gff'
include { GFF_TO_GBK                                 } from '../modules/local/gff_to_gbk'

workflow GENE_CALLING {

    take:
    assemblies_raw  // channel: full rows from fromSamplesheet [meta, assembly, gff, faa, gbk]
    assemblies      // channel: [meta, assembly]
    bakta_db        // channel: path (ignored when !params.bakta)

    main:
    ch_versions     = Channel.empty()
    annotations_fna = channel.empty()
    annotations_gbk = channel.empty()
    annotations_faa = channel.empty()
    annotations_gff = channel.empty()
    compliant_gbk   = channel.empty()
    compliant_gff   = channel.empty()

    // ------------------------------------------------------------------
    // Per-sample external gene calls from samplesheet columns
    // ------------------------------------------------------------------
    // Tuple: [meta, assembly, gene_calls_gff, gene_calls_faa, gene_calls_gbk]
    // When the gene_calls_gff column is absent, this channel is empty.
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

    // Any external sample that provides a GBK uses it directly (bypasses GFF_TO_GBK).
    ext_gbk = ensembl_gbk_branches.has_gbk
        .map    { meta, _assembly, _gff, _faa, gbk -> [meta, gbk] }
        .mix( ext_gc_branches.direct
            .filter { meta, _assembly, _gff, _faa, gbk -> gbk }
            .map    { meta, _assembly, _gff, _faa, gbk -> [meta, gbk] } )
    ext_direct_compliant_gff = ext_gc_branches.direct
        .filter { meta, _assembly, _gff, _faa, gbk -> gbk }
        .map    { meta, _assembly, gff, _faa, _gbk -> [meta, gff] }

    // ------------------------------------------------------------------
    // Kingdom lookup (needed for Prokka/Bakta routing and tRNA downstream)
    // ------------------------------------------------------------------
    LOOKUP_KINGDOM( assemblies )
    ch_versions = ch_versions.mix(LOOKUP_KINGDOM.out.versions.first())

    assemblies_with_kingdom = assemblies.join( LOOKUP_KINGDOM.out.detected_kingdom )

    // ------------------------------------------------------------------
    // Gene calling: --gene_calls (directory), or Prokka / Bakta
    // ------------------------------------------------------------------
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

        // Preserve the original samplesheet meta (needed for downstream joins)
        assemblies_by_prefix = assemblies.map { meta, assembly -> [meta.prefix, meta] }

        // Load annotation files from the appropriate tool's output directory.
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
                archaea:  it[2] == "Archaea"
            }.set { assemblies_to_annotate }

            BAKTA_BAKTA( assemblies_to_annotate.bacteria, bakta_db )
            ch_versions = ch_versions.mix(BAKTA_BAKTA.out.versions.first())

            prokka_input = assemblies_to_annotate.archaea

        } else {

            prokka_input = assemblies_for_gene_calling

        }

        // Prokka crashes with long contig names — temporarily rename before running
        SHORTEN_CONTIG_NAMES(
            prokka_input.map { meta, fasta, _kingdom -> [ meta, fasta ] },
            []
        )
        ch_versions = ch_versions.mix(SHORTEN_CONTIG_NAMES.out.versions.first())

        renamed_prokka_input = SHORTEN_CONTIG_NAMES.out.modified_fasta
            .join( prokka_input )
            .map { meta, renamed_fasta, _original_fasta, kingdom -> [ meta, renamed_fasta, kingdom ] }

        PROKKA_STANDARD( renamed_prokka_input, Channel.value("standard") )
        PROKKA_COMPLIANT( renamed_prokka_input, Channel.value("compliant") )
        ch_versions = ch_versions.mix(PROKKA_STANDARD.out.versions.first())

        // Revert contig renaming so GBK/GFF match the original FASTA
        PROKKA_STANDARD.out.gbk
            .mix(PROKKA_STANDARD.out.gff)
            .mix(PROKKA_STANDARD.out.fna)
            .groupTuple()
            .join(SHORTEN_CONTIG_NAMES.out.fasta_ids_mapping)
            .multiMap { meta, files_list, names_mapping ->
                files:   [ meta, files_list ]
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
            compliant_gbk   = compliant_gbk.mix( BAKTA_BAKTA.out.gbk ).mix( PROKKA_COMPLIANT.out.gbk )
            compliant_gff   = compliant_gff.mix( BAKTA_BAKTA.out.gff ).mix( PROKKA_COMPLIANT.out.gff )
        } else {
            annotations_gbk = REVERT_CONTIG_RENAMING.out.modified_gbk
            annotations_gff = REVERT_CONTIG_RENAMING.out.modified_gff
            annotations_fna = REVERT_CONTIG_RENAMING.out.modified_fasta
            annotations_faa = PROKKA_STANDARD.out.faa
            compliant_gbk   = PROKKA_COMPLIANT.out.gbk
            compliant_gff   = PROKKA_COMPLIANT.out.gff
        }
    }

    // Mix per-sample external gene calls into the main annotation channels.
    // Empty no-op when no external gene calls are provided.
    annotations_gff = annotations_gff.mix(ext_annotations_gff)
    annotations_faa = annotations_faa.mix(ext_annotations_faa)
    annotations_fna = annotations_fna.mix(ext_annotations_fna)
    annotations_gbk = annotations_gbk.mix(GFF_TO_GBK.out.gbk).mix(ext_gbk)
    compliant_gbk   = compliant_gbk.mix(GFF_TO_GBK.out.gbk).mix(ext_gbk)
    compliant_gff   = compliant_gff.mix(NORMALIZE_ENSEMBL_GFF.out.normalised_gff).mix(ext_direct_compliant_gff)

    emit:
    annotations_faa  = annotations_faa
    annotations_gff  = annotations_gff
    annotations_gbk  = annotations_gbk
    annotations_fna  = annotations_fna
    compliant_gbk    = compliant_gbk
    compliant_gff    = compliant_gff
    detected_kingdom = LOOKUP_KINGDOM.out.detected_kingdom
    versions         = ch_versions
}
