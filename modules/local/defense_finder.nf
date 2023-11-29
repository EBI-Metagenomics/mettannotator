process DEFENSE_FINDER {

    tag "${meta.prefix}"

    container 'biocontainers/defense-finder:1.2.0--pyhdfd78af_0'

    input:
    tuple val(meta), path(faa)
    path defense_finder_db

    output:
    tuple val(meta), path("defense_finder/${meta.prefix}_defense_finder_genes.tsv")  , emit: genes
    tuple val(meta), path("defense_finder/${meta.prefix}_defense_finder_systems.tsv"), emit: systems
    path "versions.yml"                                                              , emit: versions

    script:
    """
    defense-finder run \\
        -o defense_finder \\
        --models-dir ${defense_finder_db} \\
        ${faa}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: 1.2.0
    END_VERSIONS
    """
}