process LOOKUP_KINGDOM {

    tag "${meta.prefix}"

    label 'process_nano'

    container 'quay.io/microbiome-informatics/genomes-pipeline.python3base:v1.1'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), env("detected_kingdom"), emit: detected_kingdom
    path "versions.yml",                      emit: versions

    script:
    """
    detected_kingdom=\$(identify_kingdom.py -t ${meta.taxid} --include-kingdom)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    detected_kingdom="Bacteria"
    """
}
