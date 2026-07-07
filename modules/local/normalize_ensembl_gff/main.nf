process NORMALIZE_ENSEMBL_GFF {

    tag "${meta.prefix}"

    label 'process_nano'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(ensembl_gff)

    output:
    tuple val(meta), path("${meta.prefix}_normalised.gff"), emit: normalised_gff
    path "versions.yml",                                    emit: versions

    script:
    """
    normalize_ensembl_gff.py \\
        -i ${ensembl_gff} \\
        -o ${meta.prefix}_normalised.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.prefix}_normalised.gff
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
