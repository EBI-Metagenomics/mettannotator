/*
 * Interproscan
*/

process IPS {

    container 'quay.io/microbiome-informatics/genomes-pipeline.ips:5.62-94.0'
    containerOptions '--bind data:/opt/interproscan-5.62-94.0/data'

    label 'ips'

    input:
    tuple val(meta), path(faa_fasta)
    path interproscan_db
    path "versions.yml", emit: versions

    output:
    tuple val(meta), path('*.IPS.tsv'), emit: ips_annontations

    script:
    """
    interproscan.sh \
    -cpu ${task.cpus} \
    -dp \
    --goterms \
    -pa \
    -f TSV \
    --input ${faa_fasta} \
    -o ${meta.prefix}.IPS.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //)"
    END_VERSIONS
    """
}
