process LOOKUP_KINGDOM {

    tag "${meta.prefix}"

    label 'process_nano'

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.txt"), emit: detected_kingdom

    // For the version, I'm using the latest stable the genomes-annotation pipeline
    script:
    """
    identify_kingdom.py -t ${meta.taxid} --include-kingdom -o ${meta.prefix}_kingdom.txt
    """

    stub:
    """
    touch Bacteria.txt

    """
}
