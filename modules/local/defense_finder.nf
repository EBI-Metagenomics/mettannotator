process DEFENSE_FINDER {

    tag "${meta.prefix}"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/defense-finder:1.2.0--pyhdfd78af_0' :
        'biocontainers/defense-finder:1.2.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(faa), path(prokka_gff)
    tuple path(defense_finder_db), val(db_version)

    output:
    tuple val(meta), path("defense_finder_output/${meta.prefix}_defense_finder_genes.tsv")  , emit: genes
    tuple val(meta), path("defense_finder_output/${meta.prefix}_defense_finder_systems.tsv"), emit: systems
    tuple val(meta), path("defense_finder_output/${meta.prefix}_defense_finder.gff")        , emit: gff
    path "versions.yml"                                                                     , emit: versions

    script:
    """
    defense-finder run \\
        -o defense_finder_output \\
        --models-dir ${defense_finder_db} \\
        ${faa}

    process_defensefinder_result.py \\
        -i defense_finder_output/ \\
        -p ${prokka_gff} \\
        -o defense_finder_output/${meta.prefix}_defense_finder.gff -v 1.2.0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: 1.2.0
        defense-finder models: ${db_version}
    END_VERSIONS
    """
}
