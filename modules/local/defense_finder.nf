process DEFENSE_FINDER {

    tag "${meta.prefix}"

    // using quay.io for both singularity and docker due to and issue with the galaxy versions of the container
    // container 'biocontainers/defense-finder:2.0.0--pyhdfd78af_0'
    // container 'biocontainers/defense-finder:2.0.1--pyhdfd78af_0'
    container 'biocontainers/defense-finder:3.0.0--pyhdfd78af_0'
    //'https://depot.galaxyproject.org/singularity/defense-finder%3A3.0.0--pyhdfd78af_0'


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
    mkdir -p $PWD/home
    export HOME="$PWD/home"
    defense-finder run \\
        -o defense_finder_output \\
        --antidefensefinder \\
        --models-dir ${defense_finder_db} \\
        ${faa}

    process_defensefinder_result.py \\
        -i defense_finder_output/ \\
        -p ${prokka_gff} \\
        -o defense_finder_output/${meta.prefix}_defense_finder.gff -v 3.0.0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        defense-finder: 3.0.0
        defense-finder models: ${db_version}
    END_VERSIONS
    """
}
