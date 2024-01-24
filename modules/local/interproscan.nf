/*
 * Interproscan
*/

process INTERPROSCAN {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/genomes-pipeline.ips:5.62-94.0'
    containerOptions {
        if (workflow.containerEngine == 'singularity') {
            return "--bind ${interproscan_db}/data:/opt/interproscan-5.62-94.0/data"
        } else {
            return "-v ${interproscan_db}/data:/opt/interproscan-5.62-94.0/data"
        }
    }

    input:
    tuple val(meta), path(faa_fasta)
    tuple path(interproscan_db), val(db_version)

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
        InterProScan database: $db_version
    END_VERSIONS
    """
}
