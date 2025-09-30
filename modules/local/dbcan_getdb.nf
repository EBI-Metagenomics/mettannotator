process DBCAN_GETDB {

    tag "DBCan v5-2_9-13-2025"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "${params.dbs}", mode: 'copy'

    output:
    tuple path("dbcan/", type: "dir"), val("v5-2_9-13-2025"), emit: dbcan_db


    script:
    """
    wget -r -np -nH --cut-dirs=3 --reject "index.html*" https://bcb.unl.edu/dbCAN2/download/run_dbCAN_database_total/db_v5-2_9-13-2025/

    mv db_v5-2_9-13-2025 dbcan/

    echo 'v5-2_9-13-2025' > dbcan/VERSION.txt

    """
}
