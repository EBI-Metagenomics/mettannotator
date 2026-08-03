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
    mkdir -p dbcan

    # CAZyme databases
    wget -O dbcan/CAZy.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/CAZy.dmnd"
    wget -O dbcan/dbCAN.hmm "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/dbCAN.hmm"
    wget -O dbcan/dbCAN-sub.hmm "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/dbCAN_sub.hmm"
    wget -O dbcan/fam-substrate-mapping.tsv "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/fam-substrate-mapping.tsv"

    # CGC-related databases
    wget -O dbcan/TCDB.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/TCDB.dmnd"
    wget -O dbcan/TF.hmm "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/TF.hmm"
    wget -O dbcan/TF.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/TF.dmnd"
    wget -O dbcan/STP.hmm "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/STP.hmm"
    wget -O dbcan/PUL.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/PUL.dmnd"
    wget -O dbcan/dbCAN-PUL.xlsx "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/dbCAN-PUL.xlsx"
    wget -O dbcan/dbCAN-PUL.tar.gz "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/dbCAN-PUL.tar.gz"
    wget -O dbcan/peptidase_db.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/peptidase_db.dmnd"
    wget -O dbcan/sulfatlas_db.dmnd "https://pro.unl.edu/dbCAN2/download_file.php?file=run_dbCAN_database_total/db_v5-2_9-13-2025/sulfatlas_db.dmnd"

    tar -C dbcan -xzf dbcan/dbCAN-PUL.tar.gz


    # there is a mismatch between the expected name of the sub db and the way
    # it is named in the database the code below is a temporary fix for this
    # until it is properly addressed by the tool developers
    if [ -f "dbcan/dbCAN_sub.hmm" ]; then
        mv "dbcan/dbCAN_sub.hmm" "dbcan/dbCAN-sub.hmm"
    fi

    echo 'v5-2_9-13-2025' > dbcan/VERSION.txt

    """
}
