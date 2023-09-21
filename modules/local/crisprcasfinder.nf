process CRISPRCAS_FINDER {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.crisprcasfinder:4.3.2'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("crisprcasfinder_results/${meta.prefix}_crisprcasfinder.gff"), emit: gff
    tuple val(meta), path("crisprcasfinder_results/${meta.prefix}_crisprcasfinder.tsv"), emit: tsv
    tuple val(meta), path("crisprcasfinder_results/${meta.prefix}_crisprcasfinder_hq.gff"), emit: hq_gff

    script:
    """
    CRISPRCasFinder.pl -i $fasta \
    -so /opt/CRISPRCasFinder/sel392v2.so \
    -def G \
    -drpt /opt/CRISPRCasFinder/supplementary_files/repeatDirection.tsv \
    -outdir crisprcasfinder_results

    echo "Running post-processing"

    process_crispr_results.py \
    --tsv-report crisprcasfinder_results/TSV/Crisprs_REPORT.tsv \
    --gffs crisprcasfinder_results/GFF/*gff \
    --tsv-output crisprcasfinder_results/${meta.prefix}_crisprcasfinder.tsv \
    --gff-output crisprcasfinder_results/${meta.prefix}_crisprcasfinder.gff \
    --gff-output-hq crisprcasfinder_results/${meta.prefix}_crisprcasfinder_hq.gff \
    --fasta $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CRISPRCasFinder: \$(echo \$(CRISPRCasFinder.pl -v) | grep -o "version [0-9.]*" | sed "s/version //g")
        process_crispr_results.py: \$(process_crispr_results.py --version)
    END_VERSIONS
    """
}
