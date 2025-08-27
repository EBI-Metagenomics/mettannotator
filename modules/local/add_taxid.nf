process ADD_TAXID_TO_PROTEIN_FASTA {

    tag "${meta.prefix}"

    label 'process_nano'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("*_with_taxid.faa"), emit: annotations_faa_with_taxid
    path "versions.yml"                      , emit: versions

    script:
    """
    # Add user-provided taxid to the protein fasta file before running InterProScan on it.
    # This step is needed to be able to then run UniFIRE on the InterProScan output.

    add_taxid_to_faa.py -i ${faa} -t ${meta.taxid} -o ${meta.prefix}_with_taxid.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.prefix}_with_taxid.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
