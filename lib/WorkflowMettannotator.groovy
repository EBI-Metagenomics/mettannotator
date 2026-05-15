//
// This file holds several functions specific to the workflow/mettannotator.nf in the ebi-metagenomics/mettannotator pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowMettannotator {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, workflow) {
        if (params.add_slow_tools) {
           checkVersionConsistency(params, log, workflow.manifest.name as String, workflow.manifest.version as String)
           validateFastRunInputs(params, log)
        }
    }
    

    //
    // Verify that the fast-run results were produced by the same pipeline version.
    // by reading the version recorded in pipeline_info/software_versions.yml.
    //
    public static void checkVersionConsistency(params, log, String pipelineName, String currentVersion) {
        def versionsFile = new File("${params.add_slow_tools}/pipeline_info/software_versions.yml")
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
   
    //
    // Extract the pipeline version from the Workflow section of software_versions.yml.
    // Returns null if the key is not found.
    //
    private static String extractPipelineVersion(String yamlContent, String pipelineName) {
        def escapedName = java.util.regex.Pattern.quote(pipelineName)
        def matcher = yamlContent =~ /(?m)^\s+${escapedName}:\s*['"]?(\S+?)['"]?\s*$/
        return matcher ? (matcher[0][1] as String) : null
    }

    //
    // Required files that must exist in the fast-run output for every sample.
    //
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
    // allows the slow run to proceed even when a fast-run tool produced no output
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

    /**
     * Detect which annotator was used for a given sample by checking which
     * output subdirectory is present on disk.
     */
    private static Map<String, String> detectAnnotatorPatterns(File sampleDir) {
        def baktaDir  = new File(sampleDir, "functional_annotation/bakta")
        def prokkaDir = new File(sampleDir, "functional_annotation/prokka")

        if (baktaDir.exists() && baktaDir.isDirectory()) {
            return WorkflowMettannotator.REQUIRED_BAKTA_PATTERNS
        }
        if (prokkaDir.exists() && prokkaDir.isDirectory()) {
            return WorkflowMettannotator.REQUIRED_PROKKA_PATTERNS
        }
        return null
    }

    /**
     * Validate that all required fast-run output files are present for every
     * sample declared in the samplesheet.
     */
    public static void validateFastRunInputs(params, log) {

        def resultsDir = new File(params.add_slow_tools as String)
        def inputFile  = new File(params.input as String)

        if (!resultsDir.exists()) {
            Nextflow.error(
                "Directory specified by '--add-slow-tools' does not exist: ${resultsDir.absolutePath}"
            )
        }

        if (!resultsDir.isDirectory()) {
            Nextflow.error(
                "Path specified by '--add-slow-tools' is not a directory: ${resultsDir.absolutePath}"
            )
        }

        def lines = inputFile.readLines()
        def header = lines[0].split(',')*.trim()
        def prefixIdx = header.indexOf('prefix')

        if (prefixIdx < 0) {
            Nextflow.error("Samplesheet '${inputFile}' has no 'prefix' column.")
        }

        def errors       = []
        def warnings     = []
        def seenPrefixes = [] as Set

        lines.drop(1).eachWithIndex { line, idx ->
            if (line.trim().isEmpty()) return

            def fields    = line.split(',')*.trim()
            def prefix    = fields[prefixIdx]
            def sampleDir = new File(resultsDir, prefix)

            // Duplicate prefix guard: channel.join() on a non-unique key pairs
            // the wrong files across samples, causing silent data corruption.
            if (!seenPrefixes.add(prefix)) {
                errors << (
                    "Duplicate prefix '${prefix}' at samplesheet line ${idx + 2}.\n" +
                    "  Each sample must have a unique prefix. Duplicate prefixes cause\n" +
                    "  channel joins to silently pair the wrong files across samples."
                )
                return
            }

            if (!sampleDir.exists() || !sampleDir.isDirectory()) {
                errors << (
                    "Sample '${prefix}': output directory not found in fast-run results.\n" +
                    "  Expected: ${sampleDir.absolutePath}\n" +
                    "  Ensure the previous --fast run completed successfully and that\n" +
                    "  '--add-slow-tools' points to the correct output directory."
                )
                return
            }

            // Detect annotator and validate required files
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
 
            def annotatorName = (requiredPatterns == WorkflowMettannotator.REQUIRED_BAKTA_PATTERNS) ? "Bakta" : "Prokka"
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
                "  ${params.add_slow_tools}",
                "",
            ]

            errors.eachWithIndex { err, i ->
                msg << "  [${i + 1}] ${err.replace('\n', '\n       ')}"
                msg << ""
            }

            msg << header_line
            Nextflow.error(msg.join('\n'))
        }

        log.info "Fast-run input validation passed for all samples in '${params.add_slow_tools}'."
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
                if (!dir.exists() || !dir.isDirectory()) {
                    return []
                }

                dir.listFiles().findAll { file ->
                    file.name ==~ regex && (isLast || file.isDirectory())
                }
            }
        }

        return current
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //
    public static String toolCitationText(params) {

        def citation_text = [
                "Tools used in the workflow included:",
                "MultiQC (Ewels et al. 2016)",
                "."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        def reference_text = [
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        def methods_text = mqc_methods_yaml.text

        def engine = new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }
}