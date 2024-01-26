process LOOKUP_KINGDOM {

    tag "${meta.prefix}"

    label 'process_light'

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta)

    output:
    tuple val(meta), path("${meta_prefix}_kingdom.txt"), emit: detected_kingdom

    // For the version, I'm using the latest stable the genomes-annotation pipeline
    script:
    """
    identify_kingdom.py -t ${meta.taxid} -o ${meta.prefix}_kingdom.txt
    """

    stub:
    """
    touch ${meta.prefix}_kingdom.txt

    """
}
