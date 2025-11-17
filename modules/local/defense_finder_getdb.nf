
process DEFENSE_FINDER_GETDB {

    tag "Defense Finder Models 2.0.2"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "$params.dbs", mode: 'copy'

    output:
    tuple path("defense_finder/", type: "dir"), val("2.0.2"), emit: defense_finder_db

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/defense-finder/defense-finder-models_2.0.2.tar.gz

    tar -xvzf defense-finder-models_2.0.2.tar.gz

    mv 2.0.2/ defense_finder/

    echo '2.0.2' > defense_finder/VERSION.txt

    rm defense-finder-models_2.0.2.tar.gz
    """
}
