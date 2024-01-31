process INTEPROSCAN_GETDB {

    tag "IPRS Scan 5.65-97.0"

    container 'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9'

    publishDir "$params.dbs/", mode: 'copy'

    output:
    path "interproscan-5.65-97.0/data", type: 'dir', emit: interproscan_db

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/eggnog_db.tgz

    tar -xvzf interproscan-data-5.62-94.0.tar.gz

    rm interproscan-data-5.62-94.0.tar.gz
    """
}
