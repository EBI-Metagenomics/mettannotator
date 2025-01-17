process DBCAN_GETDB {

    tag "DBCan 4.1.3_V12"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "${params.dbs}", mode: 'copy'

    output:
    tuple path("dbcan/", type: "dir"), val("4.1.3_V12"), emit: dbcan_db


    script:
    """
    mkdir -p dbcan_db

    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/dbcan/dbcan_4.1.3_V12.tar.gz

    tar -xvzf dbcan_4.1.3_V12.tar.gz

    mv 4.1.3-V12 dbcan/

    echo '4.1.3_V12' > dbcan/VERSION.txt

    rm dbcan_4.1.3_V12.tar.gz
    """
}
