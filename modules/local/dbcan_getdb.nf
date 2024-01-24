process DBCAN_GETDB {

    tag "DBCan 4.0"

    container 'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9'

    publishDir "${params.dbs}", mode: 'copy'

    output:
    tuple path("dbcan/", type: "dir"), val("4.0"), emit: dbcan_db


    script:
    """
    mkdir -p dbcan_db

    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tools-reference-dbs/dbcan/dbcan_4.0.tar.gz

    tar -xvzf dbcan_4.0.tar.gz

    mv 4.0 dbcan/

    echo '4.0' > dbcan/VERSION.txt

    rm dbcan_4.0.tar.gz
    """
}
