process INTEPROSCAN_GETDB {

    tag "IPRS Scan 5.62-94.0"

    container 'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9'

    publishDir "${params.dbs}", mode: 'copy'

    output:
    tuple path("interproscan", type: "dir"), val("5.62-94.0"), emit: interproscan_db

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.62-94.0/interproscan-5.62-94.0-64-bit.tar.gz

    tar -xvzf interproscan-5.62-94.0-64-bit.tar.gz

    mv interproscan-5.62-94.0/ interproscan

    echo '5.62-94.0' > interproscan/VERSION.txt

    rm interproscan-5.62-94.0-64-bit.tar.gz
    """
}
