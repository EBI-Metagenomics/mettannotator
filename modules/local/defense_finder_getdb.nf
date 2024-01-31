
process DEFENSE_FINDER_GETDB {

    tag "Defense Finder Models 1.2.3"

    container 'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9'

    publishDir "$params.dbs/defense_finder_db", mode: 'copy'

    output:
    path "defense-finder-models_1.2.3/", type: "dir", emit: defense_finder_db

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tools-reference-dbs/defense-finder/defense-finder-models_1.2.3.tar.gz

    tar -xvzf defense-finder-models_1.2.3.tar.gz

    rm defense-finder-models_1.2.3.tar.gz
    """
}
