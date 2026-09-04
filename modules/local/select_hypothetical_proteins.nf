process SELECT_HYPOTHETICAL_PROTEINS {

    tag "${meta.prefix}"

    label 'process_nano'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${meta.prefix}_hypothetical.faa"), emit: hypothetical_faa

    script:
    """
    awk '/^>/ {keep = /hypothetical protein/} keep' ${faa} > ${meta.prefix}_hypothetical.faa
    """

    stub:
    """
    touch ${meta.prefix}_hypothetical.faa
    """
}
