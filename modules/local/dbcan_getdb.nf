process DBCAN_GETDB {

    tag "DBCan v5-2_9-13-2025"

    container 'docker://ubuntu:24.04'

    publishDir "${params.dbs}", mode: 'copy'

    output:
    tuple path("dbcan/"), val("v5-2_9-13-2025"), emit: dbcan_db

    shell:
    '''
    set -euo pipefail

    mkdir -p dbcan

    curl -sL \
    "https://pro.unl.edu/dbCAN2/browse_download.php?path=run_dbCAN_database_total/db_v5-2_9-13-2025" \
    | grep -o 'download_file.php?file=[^"]*' \
    | while read -r file; do
        wget -P dbcan --content-disposition "https://pro.unl.edu/dbCAN2/$file"
    done

    if [ -f "dbcan/dbCAN_sub.hmm" ]; then
        mv dbcan/dbCAN_sub.hmm dbcan/dbCAN-sub.hmm
    fi

    echo 'v5-2_9-13-2025' > dbcan/VERSION.txt
    '''
}
