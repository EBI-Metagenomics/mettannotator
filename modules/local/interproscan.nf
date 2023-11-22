/*
 * Interproscan
*/

process IPS {

    container 'quay.io/microbiome-informatics/genomes-pipeline.ips:5.62-94.0'
    containerOptions "${ workflow.containerEngine == 'singularity' ?
        '--bind data:/opt/interproscan-5.62-94.0/data':
        '-v /host/path/to/data:/opt/interproscan-5.62-94.0/data' }"

    input:
    tuple val(meta), path(faa_fasta)
    path interproscan_db

    output:
    tuple val(meta), path('*.IPS.tsv'), emit: ips_annotations
    path "versions.yml"               , emit: versions

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
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
    END_VERSIONS
    """
}
