process LOOKUP_KINGDOM {

    tag "${meta.prefix}"

    label 'process_nano'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'oras://community.wave.seqera.io/library/pip_requests_retry:4b6e4901d175ea72' :
        'community.wave.seqera.io/library/pip_requests_retry:d1a6734506332b90' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), env(detected_kingdom), emit: detected_kingdom  // noqa
    path "versions.yml",                    emit: versions

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
