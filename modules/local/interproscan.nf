/*
 * Interproscan
*/

process INTERPROSCAN {

    tag "${meta.prefix}"

    container 'quay.io/microbiome-informatics/interproscan:5.74-105.0'
    containerOptions {
        if (workflow.containerEngine in ['singularity', 'apptainer']) {
            return "--bind ${interproscan_db}/data:/opt/interproscan/data"
        } else {
            return "-v ./${interproscan_db}/data:/opt/interproscan/data"
        }
    }

    input:
    tuple val(meta), path(faa_fasta)
    tuple path(interproscan_db), val(db_version)

    output:
    tuple val(meta), path('*.IPS.tsv'), emit: ips_annotations
    tuple val(meta), path('*.IPS-raw-output.xml'), emit: ips_xml
    path "versions.yml"               , emit: versions

    script:
    """
    # Enforce encoding to prevent errors from non-ASCII characters in FASTA headers
    export JAVA_TOOL_OPTIONS="-Dfile.encoding=UTF-8"

    interproscan.sh \
    -cpu ${task.cpus} \
    -dp \
    --goterms \
    -pa \
    --input ${faa_fasta} \
    --output-file-base ${meta.prefix}.IPS-raw-output

    reformat_ips_output.py -i ${meta.prefix}.IPS-raw-output.tsv -o ${meta.prefix}.IPS.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        InterProScan: \$(interproscan.sh --version | grep -o "InterProScan version [0-9.-]*" | sed "s/InterProScan version //")
        InterProScan database: $db_version
    END_VERSIONS
    """
}
