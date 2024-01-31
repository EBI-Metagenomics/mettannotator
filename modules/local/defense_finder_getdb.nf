
process DEFENSE_FINDER_GETDB {

    tag "Defense Finder Models 1.2.3"

    container 'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9'

    publishDir "$params.dbs", mode: 'copy'

    output:
    tuple path("defense_finder/", type: "dir"), val("1.2.3"), emit: defense_finder_db

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tools-reference-dbs/defense-finder/defense-finder-models_1.2.3.tar.gz

    tar -xvzf defense-finder-models_1.2.3.tar.gz

    mv 1.2.3/ defense_finder/

    echo '1.2.3' > defense_finder/VERSION.txt

    rm defense-finder-models_1.2.3.tar.gz
    """
}
