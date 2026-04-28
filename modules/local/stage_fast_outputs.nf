/**
 * STAGE_FAST_OUTPUTS
 *
 * Stages fast-run per-sample results into the current run's output directory
 * so that slow-tool outputs land alongside them in a unified directory tree.
 *
 * Staging strategy (in priority order)
 * ──────────────────────────────────────
 * 1. cp -rl  — hard-link copy: O(1) disk space, near-instant.
 *              Works when fast-run outdir and --outdir are on the same
 *              filesystem (typical HPC/local scenario).
 * 2. rsync   — byte-copy fallback for cross-device/cloud scenarios.
 */

process STAGE_FAST_OUTPUTS {

    tag "${meta.prefix}"

    label 'process_single'

    input:
    tuple val(meta), path(fast_sample_dir)

    output:
    tuple val(meta), path("staged/**"), optional: true, emit: staged_files
 
    script:
    def prefix = meta.prefix
    """
    mkdir -p "staged"

    # Hard-link copy first (O(1) space, same filesystem)
    if cp -rl "${fast_sample_dir}/." "staged/" 2>/dev/null; then
        echo "[STAGE_FAST_OUTPUTS] Hard-linked fast-run outputs for '${prefix}'."
    else
        echo "[STAGE_FAST_OUTPUTS] Hard-link failed (cross-device?); falling back to rsync."
        rsync -a --copy-links "${fast_sample_dir}/" "staged/"
    fi

    # Drop the old merged_gff — ANNOTATE_GFF produces the definitive version
    # that integrates both fast and slow tool results.
    rm -rf "staged/functional_annotation/merged_gff"
    """

    stub:
    """
    mkdir -p "staged/functional_annotation/prokka"
    touch "staged/functional_annotation/prokka/${meta.prefix}.faa"
    """
}
