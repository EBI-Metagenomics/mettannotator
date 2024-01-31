process DBCAN_GETDB {

    tag "DBCan 4.0"

    container 'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9'

    publishDir "$params.dbs/", mode: 'copy'

    output:
    path "dbcan_4.0", type: 'dir', emit: dbcan_db

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tools-reference-dbs/dbcan/dbcan_4.0.tar.gz

    tar -xvzf dbcan_4.0.tar.gz

    rm dbcan_4.0.tar.gz
    """
}
