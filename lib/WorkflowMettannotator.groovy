//
// This file holds several functions specific to the workflow/mettannotator.nf in the ebi-metagenomics/mettannotator pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowMettannotator {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (params.add_slow_tools) {
            validateFastRunInputs(params, log)
        }
    }

    //
    // Required files that must exist in the fast-run output for every sample.
    //
    static final Map<String, String> REQUIRED_FAST_RUN_PATTERNS = [
        "Prokka protein FASTA (.faa)"  : "functional_annotation/prokka/*.faa",
        "Prokka GenBank (.gbk)"        : "functional_annotation/prokka/*.gbk",
        "Prokka GFF (.gff)"            : "functional_annotation/prokka/*.gff",
    ].asImmutable()

    //
    // Per-sample files that are expected but optional.
    //
    static final Map<String, String> EXPECTED_FAST_RUN_PATTERNS = [
        "EggNOG annotations"           : "functional_annotation/eggnog_mapper/*.emapper.annotations",
        "ncRNA table"                  : "rnas/ncrna/*.ncrna.deoverlap.tbl",
        "tRNA GFF"                     : "rnas/trna/*_trna.gff",
        "CRISPRCasFinder HQ GFF"       : "mobilome/crisprcas_finder/*_crisprcasfinder_hq.gff",
        "AMRFinder+ TSV"               : "antimicrobial_resistance/amrfinder_plus/*_amrfinderplus.tsv",
        "antiSMASH GFF"                : "biosynthetic_gene_clusters/antismash/*_antismash.gff",
        "GECCO clusters GFF"           : "biosynthetic_gene_clusters/gecco/*_gecco_clusters.gff",
        "dbCAN GFF"                    : "functional_annotation/dbcan/*_dbcan.gff",
        "DefenseFinder GFF"            : "antiphage_defense/defense_finder/*_defense_finder.gff",
        "Pseudofinder processed GFF"   : "functional_annotation/pseudofinder/*_processed_pseudogenes.gff",
    ].asImmutable()

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

        def errors   = []
        def warnings = []

        lines.drop(1).eachWithIndex { line, idx ->
            if (line.trim().isEmpty()) return

            def fields = line.split(',')*.trim()
            def prefix    = fields[prefixIdx]
            def sampleDir = new File(resultsDir, prefix)

            if (!sampleDir.exists() || !sampleDir.isDirectory()) {
                errors << (
                    "Sample '${prefix}': output directory not found in fast-run results.\n" +
                    "  Expected: ${sampleDir.absolutePath}\n" +
                    "  Ensure the previous --fast run completed successfully and that\n" +
                    "  '--add-slow-tools' points to the correct output directory."
                )
                return
            }

            REQUIRED_FAST_RUN_PATTERNS.each { label, pattern ->
                def matches = findMatches(sampleDir, pattern)

                if (matches.isEmpty()) {
                    errors << (
                        "Sample '${prefix}': required file missing — ${label}\n" +
                        "  Pattern: ${sampleDir.absolutePath}/${pattern}\n" +
                        "  The fast run may have failed or the output directory structure\n" +
                        "  does not match what mettannotator expects."
                    )
                }
            }

            EXPECTED_FAST_RUN_PATTERNS.each { label, pattern ->
                def matches = findMatches(sampleDir, pattern)

                if (matches.isEmpty()) {
                    warnings << "Sample '${prefix}': expected file not found (will be skipped) — ${label} [${sampleDir.absolutePath}/${pattern}]"
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