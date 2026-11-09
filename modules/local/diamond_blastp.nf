process DIAMOND_BLASTP {

    tag "${meta.prefix} ${db_label}"

    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/diamond:2.2.5--he361c42_0' :
        'biocontainers/diamond:2.2.5--he361c42_0' }"

    input:
    tuple val(meta), path(query_faa)
    tuple path(diamond_db), val(db_version)
    val db_label

    output:
    tuple val(meta), path("${meta.prefix}_${db_label}.tsv"), emit: hits_tsv
    path "versions.yml",                                     emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    diamond blastp \\
        --db ${diamond_db} \\
        --query ${query_faa} \\
        --out ${meta.prefix}_${db_label}.tsv \\
        --outfmt 6 qseqid stitle qlen slen length qcovhsp pident evalue bitscore \\
        --threads ${task.cpus} \\
        --header simple \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | sed 's/^.*diamond version //')
        ${db_label}: ${db_version}
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.prefix}_${db_label}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | sed 's/^.*diamond version //')
        ${db_label}: ${db_version}
    END_VERSIONS
    """
}
