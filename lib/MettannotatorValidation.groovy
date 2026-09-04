import nextflow.Nextflow

class MettannotatorValidation {

    static final Map<String, String> REQUIRED_PROKKA_PATTERNS = [
        "Prokka protein FASTA (.faa)"  : "functional_annotation/prokka/*.faa",
        "Prokka GenBank (.gbk)"        : "functional_annotation/prokka/*.gbk",
        "Prokka GFF (.gff)"            : "functional_annotation/prokka/*.gff",
    ].asImmutable()

    static final Map<String, String> REQUIRED_BAKTA_PATTERNS = [
        "Bakta protein FASTA (.faa)"  : "functional_annotation/bakta/*.faa",
        "Bakta GenBank (.gbff)"       : "functional_annotation/bakta/*.gbff",
        "Bakta GFF (.gff3)"           : "functional_annotation/bakta/*.gff3",
    ].asImmutable()

    // Annotation files expected in annotate_gff.py
    static final Map<String, String> REQUIRED_ANNOTATION_INPUTS = [
        "EggNOG annotations (.emapper.annotations)" : "functional_annotation/eggnog_mapper/*.emapper.annotations",
        "ncRNA table (.ncrna.deoverlap.tbl)"         : "rnas/ncrna/*.ncrna.deoverlap.tbl",
        "tRNA GFF (_trna.gff)"                       : "rnas/trna/*_trna.gff",
    ].asImmutable()

    // Per-tool output files from the fast run.
    // By default these are treated as required: a missing file is a hard error.
    // Pass '--allow-missing-files' to downgrade missing files to warnings, which
    // allows the full run to proceed even when a fast-run tool produced no output
    // (e.g. a genome with no CRISPR arrays, or a tool that was intentionally
    // skipped).  annotate_gff.py guards every one of these with a None-check so
    // the annotation step is safely skipped when the file is absent.
    static final Map<String, String> EXPECTED_FAST_RUN_PATTERNS = [
        "CRISPRCasFinder HQ GFF"       : "mobilome/crisprcas_finder/*_crisprcasfinder_hq.gff",
        "AMRFinder+ TSV"               : "antimicrobial_resistance/amrfinder_plus/*_amrfinderplus.tsv",
        "antiSMASH GFF"                : "biosynthetic_gene_clusters/antismash/*_antismash.gff",
        "GECCO clusters GFF"           : "biosynthetic_gene_clusters/gecco/*_gecco_clusters.gff",
        "dbCAN GFF"                    : "functional_annotation/dbcan/*_dbcan.gff",
        "DefenseFinder GFF"            : "antiphage_defense/defense_finder/*_defense_finder.gff",
        "Pseudofinder processed GFF"   : "functional_annotation/pseudofinder/*_processed_pseudogenes.gff",
    ].asImmutable()

    public static void validateSamplesheetPrefixes(params, log) {
        def inputFile = new File(params.input as String)
        def lines = inputFile.readLines()
        def header = lines[0].split(',')*.trim()
        def prefixIdx = header.indexOf('prefix')
        if (prefixIdx < 0) return

        def errors       = []
        def seenPrefixes = [] as Set

        lines.drop(1).eachWithIndex { line, idx ->
            if (line.trim().isEmpty()) return

            def fields = line.split(',')*.trim()
            def prefix = fields[prefixIdx]

            // Duplicate prefix guard: channel.join() on a non-unique key pairs
            // the wrong files across samples, causing silent data corruption.
            if (!seenPrefixes.add(prefix)) {
                errors << (
                    "Duplicate prefix '${prefix}' at samplesheet line ${idx + 2}.\n" +
                    "  Each sample must have a unique prefix. Duplicate prefixes cause\n" +
                    "  channel joins to silently pair the wrong files across samples."
                )
            }
        }

        if (errors) {
            Nextflow.error(errors.join('\n\n'))
        }
    }

    //
    // Verify that the fast-run results were produced by the same pipeline version.
    // by reading the version recorded in pipeline_info/software_versions.yml.
    //
    public static void checkVersionConsistency(params, log, String pipelineName, String currentVersion) {
        def versionsFile = new File("${params.extend_to_full}/pipeline_info/software_versions.yml")
        if (!versionsFile.exists()) {
            log.warn(
                "WARNING: '${versionsFile.absolutePath}' not found in the fast-run output.\n" +
                "  Version compatibility cannot be verified. If the fast results were generated\n" +
                "  with a different pipeline version, results may be inconsistent."
            )
            return
        }

        def fastVersion = extractPipelineVersion(versionsFile.text, pipelineName)
        if (fastVersion == null) {
            log.warn(
                "WARNING: Could not determine pipeline version from '${versionsFile.absolutePath}'.\n" +
                "  Version compatibility cannot be verified."
            )
            return
        }

        if (fastVersion == currentVersion) {
            return
        }

        def msg = (
            "Pipeline version mismatch detected.\n" +
            "  Fast results generated with : v${fastVersion}\n" +
            "  Current pipeline version    : v${currentVersion}\n" +
            "  Re-run the fast step with the current pipeline version, or use\n" +
            "  '--ignore_version_mismatch' to proceed anyway."
        )

        if (params.ignore_version_mismatch) {
            log.warn("${msg}\n  Skipping version check because '--ignore_version_mismatch' was set.")
        } else {
            Nextflow.error(msg)
        }
    }

    private static String extractPipelineVersion(String yamlContent, String pipelineName) {
        def yaml = new groovy.yaml.YamlSlurper().parseText(yamlContent)
        return yaml?.Workflow?[pipelineName] as String
    }

    /**
     * Validate that all required fast-run output files are present for every
     * sample declared in the samplesheet.
     */
    public static void validateFastRunInputs(params, log) {

        def resultsDir = new File(params.extend_to_full as String)
        def inputFile  = new File(params.input as String)

        if (!resultsDir.exists()) {
            Nextflow.error(
                "Directory specified by '--extend_to_full' does not exist: ${resultsDir.absolutePath}"
            )
        }

        if (!resultsDir.isDirectory()) {
            Nextflow.error(
                "Path specified by '--extend_to_full' is not a directory: ${resultsDir.absolutePath}"
            )
        }

        def lines = inputFile.readLines()
        def header = lines[0].split(',')*.trim()
        def prefixIdx = header.indexOf('prefix')

        if (prefixIdx < 0) {
            Nextflow.error("Samplesheet '${inputFile}' has no 'prefix' column.")
        }

        def errors   = []
        def warnings = []

        lines.drop(1).eachWithIndex { line, idx ->
            if (line.trim().isEmpty()) return

            def fields    = line.split(',')*.trim()
            def prefix    = fields[prefixIdx]
            def sampleDir = new File(resultsDir, prefix)

            if (!sampleDir.isDirectory()) {
                errors << (
                    "Sample '${prefix}': output directory not found in fast-run results.\n" +
                    "  Expected: ${sampleDir.absolutePath}\n" +
                    "  Ensure the previous --fast run completed successfully and that\n" +
                    "  '--extend_to_full' points to the correct output directory."
                )
                return
            }

            def requiredPatterns = detectAnnotatorPatterns(sampleDir)

            if (requiredPatterns == null) {
                errors << (
                    "Sample '${prefix}': no annotator output directory found.\n" +
                    "  Looked for:\n" +
                    "    ${sampleDir.absolutePath}/functional_annotation/prokka/\n" +
                    "    ${sampleDir.absolutePath}/functional_annotation/bakta/\n" +
                    "  The fast run may have failed at the annotation step."
                )
                return
            }

            def annotatorName = (requiredPatterns == MettannotatorValidation.REQUIRED_BAKTA_PATTERNS) ? "Bakta" : "Prokka"
            log.debug "Sample '${prefix}': detected annotator = ${annotatorName}"

            requiredPatterns.each { label, pattern ->
                def matches = findMatches(sampleDir, pattern)
                if (matches.isEmpty()) {
                    errors << (
                        "Sample '${prefix}': required file missing — ${label}\n" +
                        "  Pattern: ${sampleDir.absolutePath}/${pattern}\n" +
                        "  Detected annotator: ${annotatorName}\n" +
                        "  The fast run may have failed or the output directory structure\n" +
                        "  does not match what mettannotator expects."
                    )
                }
            }

            // These files are required by annotate_gff.py with no None-guard.
            REQUIRED_ANNOTATION_INPUTS.each { label, pattern ->
                def matches = findMatches(sampleDir, pattern)
                if (matches.isEmpty()) {
                    errors << (
                        "Sample '${prefix}': required annotation input missing — ${label}\n" +
                        "  Pattern: ${sampleDir.absolutePath}/${pattern}\n" +
                        "  Check that the fast run completed without errors for this sample\n" +
                        "  (look for failures in the EggNOG, ncRNA or tRNA detection steps)."
                    )
                }
            }

            EXPECTED_FAST_RUN_PATTERNS.each { label, pattern ->
                def matches = findMatches(sampleDir, pattern)

                if (matches.isEmpty()) {
                    if (params.allow_missing_files) {
                        warnings << (
                            "Sample '${prefix}': file not found (step will be skipped) — ${label}\n" +
                            "  Pattern: ${sampleDir.absolutePath}/${pattern}\n" +
                            "  Missing file accepted because '--allow-missing-files' is set."
                        )
                    } else {
                        errors << (
                            "Sample '${prefix}': required file missing — ${label}\n" +
                            "  Pattern: ${sampleDir.absolutePath}/${pattern}\n" +
                            "  All fast-run tool outputs are required by default. If this file\n" +
                            "  is intentionally absent (e.g. no CRISPR arrays detected, or the\n" +
                            "  tool was deliberately skipped), re-run with '--allow-missing-files'\n" +
                            "  to treat missing optional outputs as warnings instead of errors."
                        )
                    }
                }
            }
        }

        warnings.each { log.warn it }

        if (errors) {
            def header_line = "=" * 72
            def msg = [
                "",
                header_line,
                "  FAST-RUN INPUT VALIDATION FAILED",
                header_line,
                "  ${errors.size()} error(s) found in fast-run output directory:",
                "  ${params.extend_to_full}",
                "",
            ]

            errors.eachWithIndex { err, i ->
                msg << "  [${i + 1}] ${err.replace('\n', '\n       ')}"
                msg << ""
            }

            msg << header_line
            Nextflow.error(msg.join('\n'))
        }

        log.info "Fast-run input validation passed for all samples in '${params.extend_to_full}'."
    }

    //
    // Validate --gene_calls inputs: directory must exist and contain recognisable outputs.
    //
    public static void validateGeneCallsInputs(params, log) {
        def geneCallsDir = new File(params.gene_calls as String)
        if (!geneCallsDir.isDirectory()) {
            Nextflow.error(
                "'--gene_calls' directory not found: ${geneCallsDir.absolutePath}\n" +
                "  The path must point to an existing mettannotator output directory."
            )
        }

        def inputFile = new File(params.input as String)
        def lines = inputFile.readLines()
        def header = lines[0].split(',')*.trim()
        def prefixIdx = header.indexOf('prefix')

        if (prefixIdx < 0) {
            Nextflow.error("Samplesheet '${inputFile}' has no 'prefix' column.")
        }

        def errors = []

        lines.drop(1).each { line ->
            if (line.trim().isEmpty()) return

            def prefix = line.split(',')*.trim()[prefixIdx]
            def sampleDir = new File(geneCallsDir, prefix)

            if (!sampleDir.isDirectory()) {
                errors << (
                    "Sample '${prefix}': output directory not found in gene calls directory.\n" +
                    "  Expected: ${sampleDir.absolutePath}\n" +
                    "  Ensure '--gene_calling_only' was run for this sample and that\n" +
                    "  '--gene_calls' points to the correct output directory."
                )
                return
            }

            def hasProkka = new File("${sampleDir}/functional_annotation/prokka").isDirectory()
            def hasBakta  = new File("${sampleDir}/functional_annotation/bakta").isDirectory()

            if (!hasProkka && !hasBakta) {
                errors << (
                    "Sample '${prefix}': no gene calling output found.\n" +
                    "  Looked for:\n" +
                    "    ${sampleDir.absolutePath}/functional_annotation/prokka/\n" +
                    "    ${sampleDir.absolutePath}/functional_annotation/bakta/"
                )
            }
        }

        if (errors) {
            Nextflow.error(
                "Gene calls validation failed:\n" + errors.collect { "  - ${it}" }.join("\n")
            )
        }

        log.info "Gene calls input validation passed for all samples in '${params.gene_calls}'."
    }


    /**
     * Detect which annotator was used for a given sample by checking which
     * output subdirectory is present on disk.
     */
    private static Map<String, String> detectAnnotatorPatterns(File sampleDir) {
        def baktaDir  = new File(sampleDir, "functional_annotation/bakta")
        def prokkaDir = new File(sampleDir, "functional_annotation/prokka")

        if (baktaDir.isDirectory()) {
            return MettannotatorValidation.REQUIRED_BAKTA_PATTERNS
        }
        if (prokkaDir.isDirectory()) {
            return MettannotatorValidation.REQUIRED_PROKKA_PATTERNS
        }
        return null
    }

    /**
     * Dependency-free glob matcher for the simple relative patterns used above.
     * Supports path components separated by '/' and '*' wildcards in filenames.
     *
     * Examples supported:
     *   functional_annotation/prokka/*.faa
     *   rnas/trna/*_trna.gff
     */
    private static List<File> findMatches(File sampleDir, String pattern) {
        def parts = pattern.split('/') as List<String>
        def current = [sampleDir]

        parts.eachWithIndex { part, i ->
            def isLast = i == parts.size() - 1

            def regex = '^' + part
                .replace('.', '\\.')
                .replace('*', '.*') + '$'

            current = current.collectMany { dir ->
                if (!dir.isDirectory()) {
                    return []
                }

                dir.listFiles().findAll { file ->
                    file.name ==~ regex && (isLast || file.isDirectory())
                }
            }
        }

        return current
    }

    /**
     * Validate samplesheet rows that use per-sample external gene calls.
     * Called when the samplesheet contains a 'gene_calls_gff' column.
     *   - annotation_source must be set (prokka/bakta/ensembl) when gene_calls_gff is set
     *   - gene_calls_faa must be set when gene_calls_gff is set
     *   - both files must exist on disk
     */
    public static void validateExternalGeneCalls(params, log) {
        def inputFile = new File(params.input as String)
        def lines     = inputFile.readLines()
        def header    = lines[0].split(',')*.trim()

        def gffIdx    = header.indexOf('gene_calls_gff')

        // Column not present in samplesheet at all — nothing to validate.
        if (gffIdx < 0) return

        def prefixIdx = header.indexOf('prefix')
        def faaIdx    = header.indexOf('gene_calls_faa')
        def gbkIdx    = header.indexOf('gene_calls_gbk')
        def sourceIdx = header.indexOf('annotation_source')

        def valid_sources = ['prokka', 'bakta', 'ensembl'] as Set
        def errors        = []

        lines.drop(1).eachWithIndex { line, idx ->
            if (line.trim().isEmpty()) return

            def fields = line.split(',')*.trim()
            def rowNum = idx + 2  // 1-based, accounting for header row

            def prefix = prefixIdx >= 0 && prefixIdx < fields.size() ? fields[prefixIdx] : ''
            def gff    = gffIdx    >= 0 && gffIdx    < fields.size() ? fields[gffIdx]    : ''
            def faa    = faaIdx    >= 0 && faaIdx    < fields.size() ? fields[faaIdx]    : ''
            def gbk    = gbkIdx    >= 0 && gbkIdx    < fields.size() ? fields[gbkIdx]    : ''
            def source = sourceIdx >= 0 && sourceIdx < fields.size() ? fields[sourceIdx] : ''

            if (!gff) {
                // FAA, GBK, or annotation_source set without a GFF is likely a user mistake
                if (faa || gbk || source) {
                    errors << (
                        "Row ${rowNum} (prefix='${prefix}'): 'gene_calls_faa', 'gene_calls_gbk', or " +
                        "'annotation_source' is set but 'gene_calls_gff' is missing. " +
                        "'gene_calls_gff' is required when providing external gene call files."
                    )
                }
                return
            }

            if (!source || !valid_sources.contains(source)) {
                errors << (
                    "Row ${rowNum} (prefix='${prefix}'): 'gene_calls_gff' is set but " +
                    "'annotation_source' is missing or invalid (got '${source}').\n" +
                    "  Valid values: prokka, bakta, ensembl"
                )
            }

            if (!faa) {
                errors << (
                    "Row ${rowNum} (prefix='${prefix}'): 'gene_calls_gff' is set but " +
                    "'gene_calls_faa' is missing. Both must be provided together."
                )
            }

            if (gff && !new File(gff).exists()) {
                errors << (
                    "Row ${rowNum} (prefix='${prefix}'): gene_calls_gff file not found:\n" +
                    "  ${gff}"
                )
            }

            if (faa && !new File(faa).exists()) {
                errors << (
                    "Row ${rowNum} (prefix='${prefix}'): gene_calls_faa file not found:\n" +
                    "  ${faa}"
                )
            }

            if (gbk && !new File(gbk).exists()) {
                errors << (
                    "Row ${rowNum} (prefix='${prefix}'): gene_calls_gbk file not found:\n" +
                    "  ${gbk}"
                )
            }
        }

        if (errors) {
            Nextflow.error(
                "External gene-calls samplesheet validation failed:\n" +
                errors.collect { "  - ${it}" }.join("\n")
            )
        }
    }
}
